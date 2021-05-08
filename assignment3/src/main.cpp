#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <igl/local_basis.h>
#include <igl/grad.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h>


/*** insert any necessary libigl headers here ***/
#include <igl/boundary_loop.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/adjacency_matrix.h>
#include <igl/sum.h>
#include <igl/diag.h>
#include <igl/speye.h>
#include <igl/repdiag.h>
#include <igl/cat.h>
#include <igl/dijkstra.h>

using namespace std;
using namespace Eigen;
using Viewer = igl::opengl::glfw::Viewer;

Viewer viewer;

// vertex array, #V x3
Eigen::MatrixXd V;

// face array, #F x3
Eigen::MatrixXi F;

// UV coordinates, #V x2
Eigen::MatrixXd UV;

// Per-face color array, #F x3
Eigen::MatrixXd colors_per_face;

// Per face area, #F x1
Eigen::VectorXd areas;

// Per Face distortion value #Fx1
Eigen::VectorXd D;

bool showingUV = false;
bool freeBoundary = false;
double TextureResolution = 10;
double threshold = 0.5;
bool angle=true, area=false, length=false; 
igl::opengl::ViewerCore temp3D;
igl::opengl::ViewerCore temp2D;

// decariation 
static void computeSurfaceGradientMatrix(SparseMatrix<double> & D1, SparseMatrix<double> & D2);
static inline void SSVD2x2(const Eigen::Matrix2d& J, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V);

void Redraw()
{
	viewer.data().clear();
	areas.setZero(F.rows(), 1);
	colors_per_face.setOnes(F.rows(), 3);
	if (!showingUV)
	{
		viewer.data().set_mesh(V, F);
		viewer.data().set_face_based(false);

		if (UV.size() != 0)
		{
			viewer.data().set_uv(TextureResolution*UV);
			viewer.data().show_texture = true;

			// here we should get the distortion value
			SparseMatrix<double>  Dx, Dy;
			computeSurfaceGradientMatrix(Dx, Dy);
			Eigen::MatrixXd J;
			D.resize(F.rows());
			J.resize(F.rows(), 4);
			J.col(0) = Dx * UV.col(0);
			J.col(1) = Dy * UV.col(0);
			J.col(2) = Dx * UV.col(1);
			J.col(3) = Dy * UV.col(1);
			for (int i = 0; i < F.rows(); i++)
			{
				Eigen::Matrix2d JJ,U, S, V;
				JJ << J(i, 0), J(i, 1), J(i, 2), J(i, 3);
				SSVD2x2(JJ, U, S, V); // extract the sigmas for each face to calculate the distortion
			if(angle)	
				D[i] = pow(S(0, 0) - S(1, 1), 2); // angle perserving
			else if(length)	
				D[i] = pow(S(0, 0)-1,2) + pow(S(1, 1)-1, 2); // length perserving
			else if(area)
				D[i] = pow(S(0, 0)*S(1, 1)-1, 2); // area perserving
			}


			for (int i = 0; i < F.rows(); i++) {
				if (D[i] > threshold) 
					colors_per_face.row(i) << 1.f, 0.f, 0.f;
				else
					colors_per_face.row(i) << 1.f, 1.f, 1.f;
			}

		}
	}
	else
	{
		viewer.data().show_texture = false;
		viewer.data().set_mesh(UV, F);

	}
	viewer.data().set_colors(colors_per_face);

}

bool callback_mouse_move(Viewer &viewer, int mouse_x, int mouse_y)
{
	if (showingUV)
		viewer.mouse_mode = igl::opengl::glfw::Viewer::MouseMode::Translation;
	return false;
}

static void computeSurfaceGradientMatrix(SparseMatrix<double> & D1, SparseMatrix<double> & D2)
{
	MatrixXd F1, F2, F3;
	SparseMatrix<double> DD, Dx, Dy, Dz;

	igl::local_basis(V, F, F1, F2, F3);
	igl::grad(V, F, DD);

	Dx = DD.topLeftCorner(F.rows(), V.rows());
	Dy = DD.block(F.rows(), 0, F.rows(), V.rows());
	Dz = DD.bottomRightCorner(F.rows(), V.rows());

	D1 = F1.col(0).asDiagonal()*Dx + F1.col(1).asDiagonal()*Dy + F1.col(2).asDiagonal()*Dz;
	D2 = F2.col(0).asDiagonal()*Dx + F2.col(1).asDiagonal()*Dy + F2.col(2).asDiagonal()*Dz;
}
static inline void SSVD2x2(const Eigen::Matrix2d& J, Eigen::Matrix2d& U, Eigen::Matrix2d& S, Eigen::Matrix2d& V)
{
	double e = (J(0) + J(3))*0.5;
	double f = (J(0) - J(3))*0.5;
	double g = (J(1) + J(2))*0.5;
	double h = (J(1) - J(2))*0.5;
	double q = sqrt((e*e) + (h*h));
	double r = sqrt((f*f) + (g*g));
	double a1 = atan2(g, f);
	double a2 = atan2(h, e);
	double rho = (a2 - a1)*0.5;
	double phi = (a2 + a1)*0.5;

	S(0) = q + r;
	S(1) = 0;
	S(2) = 0;
	S(3) = q - r;

	double c = cos(phi);
	double s = sin(phi);
	U(0) = c;
	U(1) = s;
	U(2) = -s;
	U(3) = c;

	c = cos(rho);
	s = sin(rho);
	V(0) = c;
	V(1) = -s;
	V(2) = s;
	V(3) = c;
}


void ConvertConstraintsToMatrixForm(VectorXi indices, MatrixXd positions, Eigen::SparseMatrix<double> &C, VectorXd &d)
{
	C.resize(2 * indices.rows(), 2 * V.rows());
	C.setZero();
	for (int i = 0; i < indices.rows(); i++) {
		C.insert(i, indices(i)) = C.insert((indices.rows()) + i, (C.cols() / 2) + indices(i)) = 1;;
	}
	igl::cat(1, positions.col(0), positions.col(1), d);
}

void computeParameterization(int type)
{
	VectorXi fixed_UV_indices;
	MatrixXd fixed_UV_positions;

	SparseMatrix<double> A;
	VectorXd b;
	Eigen::SparseMatrix<double> C;
	VectorXd d;
	Eigen::SparseMatrix<double > finalMat;
	VectorXd bd;
	SparseMatrix<double> Dx, Dy;

	// Find the indices of the boundary vertices of the mesh and put them in fixed_UV_indices
	if (!freeBoundary)
	{
		// The boundary vertices should be fixed to positions on the unit disc. Find these position and
		// save them in the #V x 2 matrix fixed_UV_position.
		igl::boundary_loop(F, fixed_UV_indices);
		igl::map_vertices_to_circle(V, fixed_UV_indices, fixed_UV_positions);
	}
	else
	{
		// Fix two UV vertices. This should be done in an intelligent way. Hint: The two fixed vertices should be the two most distant one on the mesh.
		// save them in the #V x 2 matrix fixed_UV_position.
		int v1=0, v2=0;
		double glob_max = 0, loc_max = 0;
		std::set<int> s;	s.insert(1);
		std::vector<std::vector<int>> VV;
		igl::adjacency_list(F, VV);
		for (int j = 0; j < V.rows(); j++) {
			VectorXd mindist, del;
			mindist.resize(V.rows(), 1); del.resize(V.rows(), 1);
			igl::dijkstra(j, s, VV, mindist, del);
			int maxInd = 0;
			loc_max = 0;
			for (int i = 0; i < V.rows(); i++) {
				if (!isinf(mindist[i])) { // we dont want to count the inf points
					if (mindist[i] > loc_max) {
						maxInd = i;
						loc_max = mindist[i];
					}
				}
			}
			if (loc_max > glob_max){
				glob_max = loc_max;
				v1 = j;	v2 = maxInd;
			}
		}

		fixed_UV_indices.resize(2);
		fixed_UV_positions.resize(2, 2);

		fixed_UV_indices.coeffRef(0) = v1;
		fixed_UV_indices.coeffRef(1) = v2;

		fixed_UV_positions.coeffRef(0, 0) = 1.f;
		fixed_UV_positions.coeffRef(0, 1) = 1.f;
		fixed_UV_positions.coeffRef(1, 0) = -1.f;
		fixed_UV_positions.coeffRef(1, 1) = -1.f;
	}
	ConvertConstraintsToMatrixForm(fixed_UV_indices, fixed_UV_positions, C, d);
	b.setZero(2 * V.rows());

	// Find the linear system for the parameterization (1- Tutte, 2- Harmonic, 3- LSCM, 4- ARAP)
	// and put it in the matrix A.
	// The dimensions of A should be 2#V x 2#V.
	if (type == '1') {
		// Add your code for computing uniform Laplacian for Tutte parameterization
		// Hint: use the adjacency matrix of the mesh
		Eigen::SparseMatrix<double> adjacency;	igl::adjacency_matrix(F, adjacency);
		SparseVector<double> rowSum;	igl::sum(adjacency, 1, rowSum);// sum each row
		SparseMatrix<double> rowSumAsDiag;	igl::diag(rowSum, rowSumAsDiag);// Convert row sums into diagonal of sparse matrix
		SparseMatrix<double> L;		L = adjacency - rowSumAsDiag;// Build uniform laplacian

		SparseMatrix<double> zero; // a zeroed spare matrix
		zero.resize(V.rows(), V.rows()); zero.setZero();

		SparseMatrix<double> up, down;
		igl::cat(2, L, zero, up);	// up = [L 0]
		igl::cat(2, zero, L, down); // down=[0 L]
		igl::cat(1, up, down, A);	// A = [up;down]	
	}

	if (type == '2') {
		// Add your code for computing cotangent Laplacian for Harmonic parameterization
		// Use can use a function "cotmatrix" from libIGL, but ~~~~***READ THE DOCUMENTATION***~~~~
		SparseMatrix<double> L,zero;
		zero.resize(V.rows(), V.rows()); zero.setZero();

		igl::cotmatrix(V, F, L);

		SparseMatrix<double> up, down;
		igl::cat(2, L, zero, up);	// up = [L 0]
		igl::cat(2, zero, L, down); // down=[0 L]
		igl::cat(1, up, down, A);	// A = [up;down]
	}

	if (type == '3') {
		// Add your code for computing the system for LSCM parameterization
		// Note that the libIGL implementation is different than what taught in the tutorial! Do not rely on it!!

		// we want A=[D1 -D2;D2 D1;D2 D1;-D1 D2]
		computeSurfaceGradientMatrix(Dx, Dy);
		VectorXd ar;
		igl::doublearea(V, F, ar);
		for (int i = 0; i < ar.size(); i++) {
			ar[i] = sqrt(ar[i]);
		}
		Dx = ar.asDiagonal()*Dx; // including the face areas to the constraints
		Dy = ar.asDiagonal()*Dy;
		SparseMatrix<double> up, down, middle;

		igl::cat(2, Dx, SparseMatrix<double>(-1 * Dy), up); // up = [D1	-D2]
		igl::cat(1, igl::cat(2, Dy, Dx), igl::cat(2, Dy, Dx), middle); // middle = [D2 D1;D2 D1]
		down = -1 * up; // down = [-D1 D2]
		igl::cat(1, up, igl::cat(1, middle, down), A); // A = A=[D1 -D2;D2 D1;D2 D1;-D1 D2]

		A = A.transpose()*A;

	}

	if (type == '4') {
		// Add your code for computing ARAP system and right-hand side
		// Implement a function that computes the local step first
		// Then construct the matrix with the given rotation matrices
		return;
	}

	// Solve the linear system.
	// Construct the system as discussed in class and the assignment sheet

	// Use igl::cat to concatenate matrices
	// Use Eigen::SparseLU to solve the system. Refer to tutorial 3 for more detail

	// largeMatrix = [A C';C 0]
	SparseMatrix<double> up, down,zero;
	zero.resize(C.rows(), C.rows());	zero.setZero();

	igl::cat(2, A, SparseMatrix<double>(C.transpose()), up);	// up = [A C']
	igl::cat(2, C, zero, down);									// down = [C 0]
	igl::cat(1, up, down, finalMat);							// large_mat = [up;down]

	igl::cat(1, b, d, bd); //bd = [b;d]

	finalMat.makeCompressed();
	Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>>  solver;
	solver.analyzePattern(finalMat);
	solver.factorize(finalMat);
	// The solver will output a vector
	UV.resize(V.rows(), 2);
	Eigen::VectorXd vec = solver.solve(bd);
	UV.col(0) = vec.head(V.rows());
	UV.col(1) = vec.segment(V.rows(), V.rows());

}


bool callback_key_pressed(Viewer &viewer, unsigned char key, int modifiers) {
	switch (key) {
	case '1':
	case '2':
	case '3':
	case '4':
		computeParameterization(key);
		break;
	case '5':
		// Add your code for detecting and displaying flipped triangles in the
		// UV domain here
		break;
	case '+':
		TextureResolution /= 2;
		break;
	case '-':
		TextureResolution *= 2;
		break;
	case ' ': // space bar -  switches view between mesh and parameterization
		if (showingUV)
		{
			temp2D = viewer.core;
			viewer.core = temp3D;
			showingUV = false;
		}
		else
		{
			if (UV.rows() > 0)
			{
				temp3D = viewer.core;
				viewer.core = temp2D;
				showingUV = true;
			}
			else { std::cout << "ERROR ! No valid parameterization\n"; }
		}
		break;
	}
	Redraw();
	return true;
}

bool load_mesh(string filename)
{
	igl::read_triangle_mesh(filename, V, F);
	Redraw();
	viewer.core.align_camera_center(V);
	showingUV = false;

	return true;
}

bool callback_init(Viewer &viewer)
{
	temp3D = viewer.core;
	temp2D = viewer.core;
	temp2D.orthographic = true;

	return false;
}

int main(int argc, char *argv[]) {
	if (argc != 2) {
		cout << "Usage ex3_bin <mesh.off/obj>" << endl;
		load_mesh("../data/cathead.obj");
	}
	else
	{
		// Read points and normals
		load_mesh(argv[1]);
	}

	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();

		// Add new group
		if (ImGui::CollapsingHeader("Parmaterization", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// Expose variable directly ...
			ImGui::Checkbox("Free boundary", &freeBoundary);

			ImGui::Checkbox("Showing UV", &showingUV);
			// TODO: Add more parameters to tweak here...

			ImGui::InputDouble("Threshold for distortion", &threshold,0.2,0.5,"%.3f");

			ImGui::Text("Distortion Types:");
			if (ImGui::Checkbox("Angle perserving ", &angle)) {				
				length = 0;area = 0;
			}
			if (ImGui::Checkbox("Length perseving ", &length))
			{
				angle = 0;area = 0;
			}
			if (ImGui::Checkbox("Area perserving ", &area))
			{
				length = 0;angle = 0;
			}
		}
	};

	viewer.callback_key_pressed = callback_key_pressed;
	viewer.callback_mouse_move = callback_mouse_move;
	viewer.callback_init = callback_init;

	viewer.launch();
}
