#include <iostream>
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/writeOFF.h>

#include "gui_utils.h"

/*** insert any libigl headers here ***/
#include <igl/barycenter.h>
#include <igl/jet.h>
#include <igl/facet_components.h>
#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/facet_components.h>
#include <igl/is_edge_manifold.h>
#include <igl/is_vertex_manifold.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/boundary_loop.h>
#include <igl/loop.h>
#include <igl/edge_topology.h>
#include <igl/barycenter.h>
#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>

using namespace std;

enum MouseMode { NONE, FACE_SELECT, VERTEX_SELECT, TRANSLATE, ROTATE };
MouseMode mouse_mode = NONE;
bool doit = false;
int down_mouse_x = -1, down_mouse_y = -1;
// Vertex array, #V x3
Eigen::MatrixXd V;
// Face array, #F x3
Eigen::MatrixXi F;
// Vectors of indices for adjacency relations
std::vector<std::vector<int> > VF, VFi, VV;
// Integer vector of component IDs per face, #F x1
Eigen::VectorXi cid;
// Per-face color array, #F x3
Eigen::MatrixXd colors_per_face;

std::set<int> selected_faces;
int selected_v = -1;
void update_display(igl::opengl::glfw::Viewer& viewer);

bool callback_key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
	if (key == '1')
	{
		igl::vertex_triangle_adjacency(V, F, VF, VFi);
		// Warning: may consider using iterators instead of classic index increaminting 
		std::cout << "VF:" << std::endl;
		for (int i = 0; i < VF.size(); i++) {
			for (int j = 0; j < VF[i].size(); j++) {
				std::cout << VF[i][j] << " ";
			}
			std::cout << std::endl;
		}
	}

	if (key == '2')
	{
		//add your code for computing vertex to vertex relations here
		//store in VV
		std::cout << "VV:" << std::endl;
		igl::adjacency_list( F, VV);
		for (int i = 0; i < VV.size(); i++) {
			for (int j = 0; j < VV[i].size(); j++) {
				std::cout << VV[i][j] << " ";
			}
			std::cout << std::endl;
		}
	}

	if (key == '3')
	{
		viewer.data().clear();
		viewer.data().set_mesh(V, F);
		colors_per_face.setZero(F.rows(),3);
		//add your code for computing per-face connected components here
		//store the component labels in cid
		//compute colors for the faces based on components
		//store the colors in component_colors_per_face
		//set the viewer colors
		igl::facet_components(F, cid);
		// now we calculate the number of components and the size of each one
		int numberOfComp = cid.maxCoeff() + 1;
		vector<int> comp(numberOfComp);
		for (int i = 0; i < F.rows(); i++) {comp[cid[i]]++;}
		std::cout << "There are " << numberOfComp << " Compononet(s)" << std::endl;
		std::cout << "and their sizes (# of faces for each component) are as follows: " << std::endl;
		for (int i = 0; i < comp.size(); i++) { std::cout << comp[i] << "\t";}
		std::cout << std::endl;
		igl::jet(cid ,true , colors_per_face);
		viewer.data().set_colors(colors_per_face);
	}

	if (key == '4')
	{
		Eigen::MatrixXd Vout = V;
		Eigen::MatrixXi Fout = F;

		// Step 1: Add face midpoints vertices, and connect them with the face's points
		Eigen::MatrixXd BC;// #Fx1 of vertices for the barycenters
		igl::barycenter(V, F, BC);

		Vout.resize(V.rows() + F.rows(), 3); //update Vout to contain newly added vertices
		for (int i = 0; i < V.rows(); i++) { Vout.row(i) = V.row(i); }
		for (int i = 0; i < BC.rows(); i++) { Vout.row(V.rows() + i) = BC.row(i); }

		//for each face in F, create 3 faces in Fout (connecting the midpoints the old vertices)
		Fout.resize(3 * F.rows(), 3);
		for (int i = 0; i < F.rows(); i++)
		{
			// initilaize all 3 (new) faces to be the orignal face
			Fout.row(3 * i) = Fout.row(3 * i + 1) = Fout.row(3 * i + 2) = F.row(i);

			// then insert the face's midpoint as a vertex for those faces
			Fout(3 * i, 0) = Fout(3 * i + 1, 1) = Fout(3 * i + 2, 2) = V.rows() + i;
		}

		// Step 2 : Shift the old vertices
		igl::adjacency_list(F, VV);
		Eigen::RowVector3d sum;//we use sum as a 1x3 vector that holds the sum of multiple vectors
		Eigen::MatrixXd P = V; //P holds the shifted old vertices
		double an;
		for (int i = 0; i < V.rows(); i++)
		{
			//(re)initialize sum
			sum*=0;
			int n = VV[i].size();
			an = (4 - 2 * cos(2 * M_PI / n)) / 9;
			//sum on all neighbouring vertices
			for (int j = 0; j < VV[i].size(); j++) { sum += V.row(VV[i][j]); }
			// update the old vertices' position accordingly
			Vout.row(i) = (1 - an)*V.row(i) + (an/n)*sum;
		}

		Eigen::MatrixXi EV,FE,EF;
		int f1, f2; // the two faces surronding the edge
		int v1, v2,v3,v4; // the 4 relevant vertices relevant to flipping the edge
		igl::edge_topology(Vout, Fout, EV, FE, EF);
		for (int i = 0; i < EF.rows(); i++) {
			v1 = EV(i, 0);
			v2 = EV(i, 1);
			if (v1 < V.rows() && v2 < V.rows()) { //if both sides of the edge are on the old vertex, we need to flip it
				f1 = EF(i, 0);	f2 = EF(i, 1); // the two faces that are connected to the edge

				// this loop gets the two other relevant vertices by checking the adjacent faces for vertices we haven't "seen"
				for(int i=0;i<3;i++){
					if (Fout(f1, i) != v1 && Fout(f1, i) != v2) {
						v3 = Fout(f1, i);
					}
					if (Fout(f2, i) != v1 && Fout(f2, i) != v2) {
						v4 = Fout(f2, i);
					}
				}

				Fout.row(f1) << v1, v4, v3;
				Fout.row(f2) << v2, v3, v4;

			}
		}
		
		// Set up the viewer to display the new mesh
		V = Vout; F = Fout;

		update_display(viewer);
	}
    
    return false;
}

std::set<int> get_v_from_faces_idx(const Eigen::MatrixXi& F, std::set<int>& face_idx) {
    std::set<int> v_set;
    for (auto f: face_idx) {
        v_set.insert(F(f,0)); v_set.insert(F(f,1)); v_set.insert(F(f,2));
    }
    return v_set;
}

void extrude(igl::opengl::glfw::Viewer& viewer) {
	Eigen::MatrixXd Vout = V;
	Eigen::MatrixXi Fout = F;

	// Get selected faces
	Eigen::MatrixXi sF(selected_faces.size(), 3); int idx = 0;
	for (auto it = selected_faces.begin(); it != selected_faces.end(); it++) { sF.row(idx++) = F.row(*it); }
	if (sF.rows() == 0) return;
	// Assert selected faces are connected
	Eigen::VectorXi comp; igl::facet_components(sF, comp);
	if (comp.maxCoeff() != 0) { cout << "Error: Not a single connected component, #face_comp =  " << comp << endl; return; }

	// 1) Get the boundary vertices surrounding the selected faces
	std::vector<int> bnd_loop;
	igl::boundary_loop(sF, bnd_loop);
	//if bnd_loop is NULL aka all faces are chosen, we chose to terminate.
	if (bnd_loop.size() == 0) return;
	// 2) Duplicate boundary vertices
	Vout.resize(V.rows() + bnd_loop.size(), 3);
	for (int i = 0; i < V.rows(); i++) Vout.row(i) = V.row(i); // set vertices as old vertices
	for (int i = 0; i < bnd_loop.size(); i++) { Vout.row(V.rows() + i) = V.row(bnd_loop[i]); } // create new vertices as duplicates of the faces boundary

	// 3) Compute direction T: The average of selected face normals
	Eigen::RowVector3d T; T.setZero(); Eigen::MatrixXd FN;
	igl::per_face_normals(V, F, FN);
	for (auto it = selected_faces.begin(); it != selected_faces.end(); it++) { T += FN.row(*it); }
	T.normalize(); T *= 0.25*(V.row(bnd_loop[1]) - V.row(bnd_loop[0])).norm();

	// 4) Offset old vertices by T
	std::set<int> inner_v = get_v_from_faces_idx(F, selected_faces);;
	for (auto v : inner_v) {
		Vout.row(v) += T;
	}

	// 5) Update Fout 
	// Add your code for updating Fout here
	Fout.resize(F.rows() + 2 * bnd_loop.size(), 3); // 2 new faces per new edge (= per new vertex)
	for (int i = 0; i < F.rows(); i++) Fout.row(i) = F.row(i); // set first 'F.rows()' faces as the old faces

	// 5.1) Get the set of faces containing the old boundary vertices (hint: call igl::vertex_triangle_adjacency on the old 'F')
	std::set<int> Fa; // The faces that contain BOTH the in_bound and loop_bound vertices
	igl::vertex_triangle_adjacency(V, F, VF, VFi);
	for (auto v : inner_v) {
		// Warning: try something better than int i to iterate
		for (int i = 0; i < VF[v].size(); i++) {
			Fa.insert(VF[v][i]);
		}
	}

	// 5.2) Get the "outer" set of faces containing the boundary vertices 
	//      (hint: call std::set_difference to compute the difference between the previously computed set of faces, and the selected faces)
	std::set<int> Fo; // The faces that only contain the loop_bound vertices and is outside of the selected area
	std::set_difference(Fa.begin(), Fa.end(), selected_faces.begin(), selected_faces.end(), std::inserter(Fo, Fo.end()));
	
	// 5.3) Edit old outer faces indices, replacing the old vertices with the indices of the duplicated boundary vertices
	std::vector<int>::iterator res;
	for (auto ind : Fo) {
		for (int i = 0; i < 3; i++) {
			res = std::find(bnd_loop.begin(), bnd_loop.end(), Fout(ind,i));
			if (res != bnd_loop.end()) {
				Fout(ind,i) = V.rows() +  int(std::distance(bnd_loop.begin(), res));
			}
			
		}
	}
		// 5.4) Add new faces, 2 per edge
		int f_idx = F.rows();
		for (int i = 0; i < bnd_loop.size(); i++) {
			int v1, v2, v3, v4;
			v1 = bnd_loop[i];
			v2 = bnd_loop[(i + 1) % (bnd_loop.size()) ];
			v3 = (V.rows() + i);
			v4 = (V.rows() + (i + 1) % (bnd_loop.size()));
			// set v1,v2,v3,v4 correctly
			Fout.row(f_idx++) << v1, v4, v2;
			Fout.row(f_idx++) << v1, v3, v4;
		}
		
		// 6) Check that the new mesh is a manifold (call is_edge_manifold, is_vertex_manifold on Vout,Fout)
		if (!igl::is_vertex_manifold(Fout, Eigen::MatrixXi())){
			std::cout << "ERROR extruding: result is non vertex manifold" << std::endl;
			return;
		}
		if (!igl::is_edge_manifold(Fout)) {
			std::cout << "ERROR extruding: result is non edge manifold" << std::endl;
			return;
		}
		// 7) Update V,F
		V = Vout; 
		F = Fout; 

		// Update gui and move to edit-translate mode
		colors_per_face = Eigen::MatrixXd::Ones(F.rows(), 3); // number of faces has changed
		viewer.data().clear();
		viewer.data().set_mesh(V, F);
		for (auto f : selected_faces) { colors_per_face.row(f) << 1, 0, 0; }
		viewer.data().set_colors(colors_per_face);
		mouse_mode = TRANSLATE;
	}

void clear_selection(igl::opengl::glfw::Viewer& viewer) {
    selected_faces.clear();
    selected_v = -1;
    colors_per_face = Eigen::MatrixXd::Ones(F.rows(),3);
    viewer.data().clear();
    viewer.data().set_mesh(V,F);
    viewer.data().set_colors(colors_per_face);
}

void export_mesh() {
    std::string f = igl::file_dialog_save();
    igl::writeOFF(f,V,F);
}

bool callback_init(igl::opengl::glfw::Viewer& viewer)
{
    return false;
}

bool callback_mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier) 
{
    down_mouse_x = viewer.current_mouse_x;
    down_mouse_y = viewer.current_mouse_y;

    if (mouse_mode == FACE_SELECT) 
	{
        int f = pick_face(viewer, down_mouse_x, down_mouse_y,V,F);
        if (f !=-1)  
		{
            selected_faces.insert(f);
            selected_v = -1;
            // update face colors
            //colors_per_face.setConstant(colors_per_face.rows(), colors_per_face.cols(), 0);
            colors_per_face.row(f) << 1,0,0;
            viewer.data().set_colors(colors_per_face);
        }
        
    } 
	else if (mouse_mode == VERTEX_SELECT) 
	{
        int v = pick_vertex(viewer, down_mouse_x, down_mouse_y,V,F);
        if (v !=-1) 
		{
            selected_v = v;
            selected_faces.clear(); 
            update_display(viewer);
            viewer.data().set_points(V.row(selected_v),Eigen::RowVector3d(1,0,0));
        }
    } 
	else if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE)) 
	{
        if (!selected_faces.empty()) 
		{
            int f = pick_face(viewer, down_mouse_x, down_mouse_y,V,F);
            if (std::find(selected_faces.begin(),selected_faces.end(),f)!= selected_faces.end()) 
			{
                doit = true;
            }
        } 
		else if (selected_v != -1) 
		{
            int v = pick_vertex(viewer, down_mouse_x, down_mouse_y,V,F);
            if (v == selected_v) 
			{
                viewer.data().set_points(V.row(selected_v),Eigen::RowVector3d(1,0,0));
                doit = true;
            }
        }
    }
    return false;
}

Eigen::RowVector3d get_face_avg(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, std::set<int>& selected_faces) {
    
    Eigen::RowVector3d avg; avg << 0,0,0;
    std::set<int> v_set = get_v_from_faces_idx(F, selected_faces);
    for (auto v: v_set) {
        avg += V.row(v);
    }
    avg/= v_set.size();
    
    return avg;
}

bool callback_mouse_move(igl::opengl::glfw::Viewer& viewer, int mouse_x, int mouse_y) 
{
	if (!doit)
	{
		return false;
	}
    if ((mouse_mode == TRANSLATE) || (mouse_mode == ROTATE))
	{
        if (!selected_faces.empty())
		{
            Eigen::RowVector3d face_avg_pt = get_face_avg(V,F,selected_faces);
            std::set<int> v_idx = get_v_from_faces_idx(F,selected_faces);
            if (mouse_mode == TRANSLATE)
			{
                Eigen::Vector3f translation = computeTranslation(viewer,
                                             mouse_x,
                                             down_mouse_x,
                                             mouse_y,
                                             down_mouse_y,
                                             face_avg_pt);

                for (auto v_i : v_idx) 
				{
					V.row(v_i) += translation.cast<double>();
				}
            } 
			else 
			{ // ROTATE
                Eigen::Vector4f rotation = computeRotation(viewer,
                                 mouse_x,
                                 down_mouse_x,
                                 mouse_y,
                                 down_mouse_y,
                                 face_avg_pt);

                for (auto v_i : v_idx) 
				{
					Eigen::RowVector3f goalPosition = V.row(v_i).cast<float>();
					goalPosition -= face_avg_pt.cast<float>();
					Eigen::RowVector4f goalPosition4 = Eigen::RowVector4f(goalPosition[0], goalPosition[1], goalPosition[2], 0);
					igl::rotate_by_quat(goalPosition4.data(), rotation.data(), goalPosition4.data());
					goalPosition = Eigen::RowVector3f(goalPosition4[0], goalPosition4[1], goalPosition4[2]);
					goalPosition += face_avg_pt.cast<float>();
					V.row(v_i) = goalPosition.cast<double>();
                }
            }

            viewer.data().set_mesh(V,F);
            down_mouse_x = mouse_x;
            down_mouse_y = mouse_y;
            return true;    
        } 
		else if ((selected_v!=-1) && (mouse_mode == TRANSLATE)) 
		{
            Eigen::Vector3f translation = computeTranslation(viewer,
                                             mouse_x,
                                             down_mouse_x,
                                             mouse_y,
                                             down_mouse_y,
                                             V.row(selected_v));

            V.row(selected_v) += translation.cast<double>();
            viewer.data().set_mesh(V,F);
            viewer.data().set_points(V.row(selected_v),Eigen::RowVector3d(1,0,0));
            down_mouse_x = mouse_x;
            down_mouse_y = mouse_y;
            return true;
        }
    }

    return false;
}

bool callback_mouse_up(igl::opengl::glfw::Viewer& viewer, int button, int modifier) {
    doit = false;
    return false;
}

void update_display(igl::opengl::glfw::Viewer& viewer) {
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    colors_per_face = Eigen::MatrixXd::Ones(F.rows(),3);
    viewer.data().set_colors(colors_per_face);
}

int main(int argc, char *argv[]) {
	// Init the viewer
	igl::opengl::glfw::Viewer viewer;

	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	// Add content to the default menu window
	menu.callback_draw_viewer_menu = [&]()
	{
		menu.draw_viewer_menu();

		ImGui::Combo("Mouse Mode", (int *)(&mouse_mode), "None\0Face Selection\0Vertex Selection\0Translate\0Roate\0\0");

		if (ImGui::Button("Extrude - W.W was Here"))
		{
			extrude(viewer);
		}

		if (ImGui::Button("Clear Selection"))
		{
			clear_selection(viewer);
		}

		if (ImGui::Button("Export Mesh"))
		{
			export_mesh();
		}
	};

    //Viewer viewer;
    viewer.callback_key_down = callback_key_down;
    
    if (argc == 2)
    {
      // Read mesh
      igl::readOFF(argv[1],V,F);
      
    }
    else
    {
      // Read mesh
      igl::readOFF("../data/cube.off",V,F);
    }

    viewer.data().set_mesh(V,F);
    viewer.data().compute_normals();
    viewer.core.align_camera_center(V,F);
    viewer.callback_init = callback_init;
    viewer.callback_mouse_down = callback_mouse_down;
    viewer.callback_mouse_up = callback_mouse_up;
    viewer.callback_mouse_move = &callback_mouse_move;
    update_display(viewer);
    viewer.launch();
}
