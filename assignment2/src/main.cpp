#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
/*** insert any necessary libigl headers here ***/
#include <igl/per_face_normals.h>
#include <igl/copyleft/marching_cubes.h>

using namespace std;
using Viewer = igl::opengl::glfw::Viewer;

// Input: imported points, #P x3
Eigen::MatrixXd P;

// Input: imported normals, #P x3
Eigen::MatrixXd N;

// Intermediate result: constrained points, #C x3
Eigen::MatrixXd constrained_points;

// Intermediate result: implicit function values at constrained points, #C x1
Eigen::VectorXd constrained_values;

// Intermediate result: constrained point colors, for display, #C x3
Eigen::MatrixXd constrained_colors;

// Parameter: degree of the polynomial
int polyDegree = 0;

// Parameter: Wendland weight function radius (make this relative to the size of the mesh)
float wendlandRadius = 0.1;

// Parameter: grid resolution
int resolution = 20;

// Intermediate result: grid points, at which the imlicit function will be evaluated, #G x3
Eigen::MatrixXd grid_points;

// Intermediate result: implicit function values at the grid points, #G x1
Eigen::VectorXd grid_values;

// Intermediate result: grid point colors, for display, #G x3
Eigen::MatrixXd grid_colors;

// Intermediate result: grid lines, for display, #L x6 (each row contains
// starting and ending point of line segment)
Eigen::MatrixXd grid_lines;

// Output: vertex array, #V x3
Eigen::MatrixXd V;

// Output: face array, #F x3
Eigen::MatrixXi F;

// Output: face normals of the reconstructed mesh, #F x3
Eigen::MatrixXd FN;

// Spatial index queue, resolution being 100
const int reso = 1;
std::vector<int> spatialIndP[reso][reso][reso];
std::vector<int> spatialIndC[reso][reso][reso];

// Functions
void createGrid();
void evaluateImplicitFunc();
void getLines();
bool callback_key_down(Viewer& viewer, unsigned char key, int modifiers);

// My additional functions
void initQueueP();
Eigen::RowVector3d getClosest(Eigen::RowVector3d q);
std::vector<int> getNeighbours(Eigen::RowVector3d q);
void adjustAxis();

void initQueueC()
{
	Eigen::RowVector3d bb_min, bb_max;
	bb_min = P.colwise().minCoeff();
	bb_max = P.colwise().maxCoeff();
	Eigen::RowVector3d dim = bb_max - bb_min;

	Eigen::RowVector3d scalar = ((reso - 1)*dim.inverse());
	Eigen::RowVector3d Pi;
	// insert only non zero weighted constraints
	for (int i = constrained_points.rows() / 3; i < constrained_points.rows(); i++) {
		Pi = constrained_points.row(i) - bb_min;
		Pi = Pi.cwiseProduct(scalar);
		if (constrained_values(i) != 0) {
			spatialIndC[int(Pi(0))][int(Pi(1))][int(Pi(2))].push_back(i);
		}


	}
}

void initQueueP() {
	Eigen::RowVector3d bb_min, bb_max;
	bb_min = P.colwise().minCoeff();
	bb_max = P.colwise().maxCoeff();
	Eigen::RowVector3d dim = bb_max - bb_min;

	Eigen::RowVector3d scalar((reso - 1)*dim.inverse());
	for (int i = 0; i < P.rows(); i++) {
		Eigen::RowVector3d Pi = P.row(i) - bb_min;
		Pi = Pi.cwiseProduct(scalar);
		spatialIndP[int(Pi(0))][int(Pi(1))][int(Pi(2))].push_back(i);
	}

}

// Return the closest point in P to given q
Eigen::RowVector3d getClosest(Eigen::RowVector3d q) {

	// the naive solution
	// init the values well
	//Eigen::RowVector3d closest;
	//double minDist;
	//int k = 0;
	//int minInd = 0;
	//closest = P.row(0);
	//minDist = (q - closest).norm();
	//// now make sure that toAdd is the closest point to Pi
	////cont = false;
	//for (k = 0; k < P.rows(); k++) {
	//	if ((q - P.row(k)).norm() < minDist) {
	//		closest = P.row(k);
	//		minDist = (q - P.row(k)).norm();
	//		minInd = k;
	//	}
	//}
	// The following implementation uses the spatial index accelarted data structure,
	// [For time Optimiazation] DONT FORGET TO CALL createQueue()

	Eigen::RowVector3d closest;
	double minDist;
	int minInd = 0;
	minDist = 1000;
	Eigen::RowVector3d bb_min, bb_max;
	bb_min = P.colwise().minCoeff();
	bb_max = P.colwise().maxCoeff();
	Eigen::RowVector3d dim = bb_max - bb_min;
	Eigen::RowVector3d scalar((reso - 1)*dim.inverse());
	std::vector<int> curVec;
	Eigen::RowVector3d tempo = q - bb_min;
	Eigen::RowVector3i ind(int(tempo(0)*scalar(0)), int(tempo(1)*scalar(1)), int(tempo(2)*scalar(2)));
	for (int a = ind(0) - 1; a <= ind(0) + 1; a++) {
		if (a<0 || a>reso - 1) continue;
		for (int b = ind(1) - 1; b <= ind(1) + 1; b++) {
			if (b<0 || b>reso - 1) continue;
			for (int c = ind(2) - 1; c <= ind(2) + 1; c++) {
				if (c<0 || c>reso - 1) continue;

				curVec = spatialIndP[a][b][c];
				for (int i = 0; i < curVec.size(); i++) {
					if ((q - P.row(curVec[i])).norm() <= minDist) {
						closest = P.row(curVec[i]);
						minDist = (q - P.row(curVec[i])).norm();
					}
				}
			}
		}
	}

	return closest;


}

// Return a vector of indices for the Points in C within a radius of wendlandRadius from q
// *AND HAVE A NONZERO CONSTRAINED VALUE*
std::vector<int> getNeighbours(Eigen::RowVector3d q) {
	// the naive solution
	//std::vector<int> ind;
	//// we know for sure that all the 3 entries at the beggening of C
	//// have a constrained value of 0, so we skip
	//for (int i = constrained_points.rows() / 3; i < constrained_points.rows(); i++) {
	//	if ((constrained_points.row(i) - q).norm() <= wendlandRadius) {
	//		ind.push_back(i);
	//	}
	//}

	std::vector<int> curVec;
	std::vector<int> ind;
	Eigen::RowVector3d bb_min, bb_max;
	bb_min = P.colwise().minCoeff();
	bb_max = P.colwise().maxCoeff();
	Eigen::RowVector3d dim = bb_max - bb_min;
	Eigen::RowVector3d scalar((reso - 1)*dim.inverse());
	std::vector<int> point_segments(3);
	point_segments[0] = int((q(0) - bb_min(0))*scalar(0));
	point_segments[1] = int((q(1) - bb_min(1))*scalar(1));
	point_segments[2] = int((q(2) - bb_min(2))*scalar(2));

	for (int a = (point_segments[0] - int(wendlandRadius*scalar(0))); a <= (point_segments[0] + int(wendlandRadius*scalar(0))); a++) {
		if (a<0 || a>reso - 1) continue;
		for (int b = (point_segments[1] - int(wendlandRadius*scalar(1))); b <= (point_segments[1] + int(wendlandRadius*scalar(1))); b++) {
			if (b<0 || b>reso - 1) continue;
			for (int c = (point_segments[2] - int(wendlandRadius*scalar(2))); c <= (point_segments[2] + int(wendlandRadius*scalar(2))); c++) {
				if (c<0 || c>reso - 1) continue;

				curVec = spatialIndC[a][b][c];

				for (int i = 0; i < curVec.size(); i++) {
					if (((constrained_points.row(curVec[i])) - q).norm() <= wendlandRadius){
						ind.push_back(curVec[i]);
					}
				}

			}
		}
	}
	return ind;
}


// alligns the object to be centered at {0,0,0} and alligned with axis
void adjustAxis()
{
	Eigen::RowVector3d m; // the center of the points
	Eigen::Matrix3d scatter; // the scatter matrix
	scatter.resize(3, 3);

	m = P.colwise().sum() / P.rows();

	for (int i = 0; i < P.rows(); i++) {
		P.row(i) = P.row(i) - m;
	}

	scatter = P.transpose()*P;
	Eigen::EigenSolver<Eigen::MatrixXd> es(scatter);
	std::cout << "Scatter matrix: " << std::endl << scatter << std::endl;
	Eigen::MatrixXd eigenvalues = scatter.eigenvalues().real();
	eigenvalues.normalize();

	std::cout << "Eigen values " << std::endl << eigenvalues << std::endl;

	double avgEigen = eigenvalues.mean();

	eigenvalues(0) = eigenvalues(0) - avgEigen; 	eigenvalues(1) = eigenvalues(1) - avgEigen;
	eigenvalues(2) = eigenvalues(2) - avgEigen;
	eigenvalues = eigenvalues.cwiseProduct(eigenvalues);
	std::cout << "the Variances are: " << std::endl << eigenvalues << std::endl;
	if(eigenvalues(0) < 0.3 && eigenvalues(1) < 0.3 && eigenvalues(2) < 0.3){
		std::cout << "There is NO preferable direction " << std::endl;
		return;
	}

	std::cout << "Eigen vecs " << std::endl << es.eigenvectors().real() << std::endl;

	// now try and change the basis
	for (int i = 0; i < P.rows(); i++) {
		P.row(i) = es.eigenvectors().real().colPivHouseholderQr().solve(P.row(i).transpose()).transpose();
	}

}

// function to calculate wendland distance
double WendlandDist(double r) {
	return pow(1 - (r / wendlandRadius), 4)*((4 * r / wendlandRadius) + 1);
}

// Creates a grid_points array for the simple sphere example. The points are
// stacked into a single matrix, ordered first in the x, then in the y and
// then in the z direction. If you find it necessary, replace this with your own
// function for creating the grid.
void createGrid() {
    grid_points.resize(0, 3);
    grid_colors.resize(0, 3);
    grid_lines. resize(0, 6);
    grid_values.resize(0);
    V. resize(0, 3);
    F. resize(0, 3);
    FN.resize(0, 3);

    // Grid bounds: axis-aligned bounding box
    Eigen::RowVector3d bb_min, bb_max;
    bb_min = P.colwise().minCoeff();
    bb_max = P.colwise().maxCoeff();

    // Bounding box dimensions
    Eigen::RowVector3d dim = bb_max - bb_min;

    // Grid spacing
    const double dx = dim[0] / (double)(resolution - 1);
    const double dy = dim[1] / (double)(resolution - 1);
    const double dz = dim[2] / (double)(resolution - 1);
    // 3D positions of the grid points -- see slides or marching_cubes.h for ordering
    grid_points.resize(resolution * resolution * resolution, 3);
    // Create each gridpoint
    for (unsigned int x = 0; x < resolution; ++x) {
        for (unsigned int y = 0; y < resolution; ++y) {
            for (unsigned int z = 0; z < resolution; ++z) {
                // Linear index of the point at (x,y,z)
                int index = x + resolution * (y + resolution * z);
                // 3D point at (x,y,z)
                grid_points.row(index) = bb_min + Eigen::RowVector3d(x * dx, y * dy, z * dz);
            }
        }
    }
}

// Function for explicitly evaluating the implicit function for a sphere of
// radius r centered at c : f(p) = ||p-c|| - r, where p = (x,y,z).
// This will NOT produce valid results for any mesh other than the given
// sphere.
// Replace this with your own function for evaluating the implicit function
// values at the grid points using MLS
void evaluateImplicitFunc() {
    // Sphere center
    auto bb_min = grid_points.colwise().minCoeff().eval();
    auto bb_max = grid_points.colwise().maxCoeff().eval();
    Eigen::RowVector3d center = 0.5 * (bb_min + bb_max);

    double radius = 0.5 * (bb_max - bb_min).minCoeff();
    // Scalar values of the grid points (the implicit function values)
    grid_values.resize(resolution * resolution * resolution);
	int n=0; 
	std::vector<int> ind;
	Eigen::MatrixXd bs,weights;
	Eigen::VectorXd fi, coef,bC;
    // Evaluate sphere's signed distance function at each gridpoint.
    for (unsigned int x = 0; x < resolution; ++x) {
        for (unsigned int y = 0; y < resolution; ++y) {
            for (unsigned int z = 0; z < resolution; ++z) {
				// Linear index of the point at (x,y,z)
				int index = x + resolution * (y + resolution * z);
				
				// get the indices of the neighbouring non zero weight points 
				ind = getNeighbours(grid_points.row(index));
				n = ind.size();
				// if there are no points satisfying the condition, we assign a large positive number
				if (n==0) {
					grid_values[index] = INFINITE;
					continue;
				}
				// Now build the necessary matrices
				weights.resize(n, n);
				weights.setZero();
				for (int i = 0; i < n; i++) {
					weights(i,i) = WendlandDist((grid_points.row(index) - constrained_points.row(ind[i])).norm());
				}

				fi.resize(n, 1);
				for (int i = 0; i < n; i++) {
					fi[i] = constrained_values[ind[i]];
				}

				//shortcuts
				double a = grid_points(index, 0), b = grid_points(index, 1), c = grid_points(index, 2);
				// filling the coeffecients and exponents according to the polyDegree
				switch (polyDegree) {
				case 0:
					bs.resize(n, 1);
					bC.resize(1, 1);
					coef.resize(1, 1);
					bC << 1;
					for (int i = 0; i < n; i++) {
						bs.row(i) << 1;
					}
					break;
				case 1:
					bs.resize(n, 4); //4 is the number of coefecients
					bC.resize(4, 1);
					coef.resize(4, 1);
					bC << 1, a, b, c;
					for (int i = 0; i < n; i++) {
						bs.row(i) << 1, constrained_points.row(ind[i]);
					}
					break;
				case 2:
					bs.resize(n, 10); //10 is the number of coefecients
					bC.resize(10, 1);
					coef.resize(10, 1);
					// todo: validate the ordering
					for (int i = 0; i < n; i++) {
						double x = constrained_points(ind[i], 0), y = constrained_points(ind[i], 1), z = constrained_points(ind[i], 2);
						bs.row(i) << 1, x, y, z, x*x, x*y, y*y, y*z, x*z, z*z;
					}
					bC << 1, a,b,c,a*a, a*b, b*b, b*c, a*c, c*c;
					break;
				default:
					std::cout << "Unvalid value for polyDegree - must be 0,1,2 but got " << polyDegree << std::endl;
					return;
				}

				
				coef = (weights*bs).colPivHouseholderQr().solve(weights*fi);
					

				grid_values[index] = coef.dot(bC);
				ind.clear();
			}
			/*
            // Value at (x,y,z) = implicit function for the sphere
            grid_values[index] = (grid_points.row(index) - center).norm() - radius;*/
            
        }
    }
}

// Code to display the grid lines given a grid structure of the given form.
// Assumes grid_points have been correctly assigned
// Replace with your own code for displaying lines if need be.
void getLines() {
    int nnodes = grid_points.rows();
    grid_lines.resize(3 * nnodes, 6);
    int numLines = 0;

    for (unsigned int x = 0; x<resolution; ++x) {
        for (unsigned int y = 0; y < resolution; ++y) {
            for (unsigned int z = 0; z < resolution; ++z) {
                int index = x + resolution * (y + resolution * z);
                if (x < resolution - 1) {
                    int index1 = (x + 1) + y * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (y < resolution - 1) {
                    int index1 = x + (y + 1) * resolution + z * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
                if (z < resolution - 1) {
                    int index1 = x + y * resolution + (z + 1) * resolution * resolution;
                    grid_lines.row(numLines++) << grid_points.row(index), grid_points.row(index1);
                }
            }
        }
    }

    grid_lines.conservativeResize(numLines, Eigen::NoChange);
}

bool callback_key_down(Viewer &viewer, unsigned char key, int modifiers) {
    if (key == '1') {
        // Show imported points
        viewer.data().clear();
        viewer.core.align_camera_center(P);
        viewer.data().point_size = 11;
        viewer.data().add_points(P, Eigen::RowVector3d(0,0,0));
    }

    if (key == '2') {
		// Show all constraints
        viewer.data().clear();
        viewer.core.align_camera_center(P);

		initQueueP();
        // Add your code for computing auxiliary constraint points here
		int n = P.rows();

		constrained_points.resize(3 * P.rows(), 3);
		constrained_values.resize(3 * P.rows());
		// 1) For each point pi in the point cloud, add a constraint of the form f(pi) = 0.
		for (int i = 0; i < P.rows(); i++) { 
			constrained_points.row(i) = P.row(i); 
			constrained_values[i] = 0;
		}

		// 2) Fix an epsilon value, for instance epsilon = 0.01 x bounding box diagonal

		//this block is taken driectly from createGrid() method
		Eigen::RowVector3d bb_min, bb_max;
		bb_min = P.colwise().minCoeff();
		bb_max = P.colwise().maxCoeff();
		Eigen::RowVector3d dim = bb_max - bb_min;

		double EPS = 0.01*dim.sum();
		double e;
		bool cont; // boolean indicator for if toAdd is the closest point to Pi
		Eigen::RowVector3d toAdd,Pi,ni;
		for (int i = 0; i < P.rows(); i++) {
			// the constrained process for Pi+eN
			e = 2*EPS; // set up variables
			Pi = P.row(i);
			N.row(i).normalize(); ni = N.row(i); 
			do {
				e = e / 2;
				toAdd = Pi + e*ni;
			}while((getClosest(toAdd)-Pi).norm()!=0);
			constrained_points.row(i+n) = toAdd;
			constrained_values[i + n] = e;

			// the constrained process for Pi-eN
			e = 2 * EPS;
			do {
				e = e / 2;
				toAdd = Pi - e*ni;
			} while ((getClosest(toAdd) - Pi).norm() != 0);
			constrained_points.row(i+2*n) = toAdd;
			constrained_values[i + 2*n] = -e;
		}

        // Add code for displaying all points, as above
		viewer.data().clear();
		viewer.data().point_size = 6;
		constrained_colors.resize(3 * P.rows(), 3);
		for(int i=0;i<P.rows();i++){
			constrained_colors.row(i) = Eigen::RowVector3d(0, 0, 1);
			constrained_colors.row(i + n) = Eigen::RowVector3d(1, 0, 0);
			constrained_colors.row(i + 2*n) = Eigen::RowVector3d(0, 1, 0);
		}
		viewer.data().add_points(constrained_points, constrained_colors);

    }

    if (key == '3') {
        // Show grid points with colored nodes and connected with lines
        viewer.data().clear();
        viewer.core.align_camera_center(P);

		initQueueC();
        // Add code for creating a grid
        // Add your code for evaluating the implicit function at the grid points
        // Add code for displaying points and lines
        // You can use the following example:

        /*** begin: sphere example, replace (at least partially) with your code ***/
        // Make grid
        createGrid();

        // Evaluate implicit function
        evaluateImplicitFunc();

        // get grid lines
        getLines();

        // Code for coloring and displaying the grid points and lines
        // Assumes that grid_values and grid_points have been correctly assigned.
        grid_colors.setZero(grid_points.rows(), 3);

        // Build color map
        for (int i = 0; i < grid_points.rows(); ++i) {
            double value = grid_values(i);
            if (value < 0) {
                grid_colors(i, 1) = 1;
            }
            else {
                if (value > 0)
                    grid_colors(i, 0) = 1;
            }
        }
        // Draw lines and points
        viewer.data().point_size = 8;
        viewer.data().add_points(grid_points, grid_colors);
        viewer.data().add_edges(grid_lines.block(0, 0, grid_lines.rows(), 3),
                              grid_lines.block(0, 3, grid_lines.rows(), 3),
                              Eigen::RowVector3d(0.8, 0.8, 0.8));
        /*** end: sphere example ***/
    }

    if (key == '4') {
        // Show reconstructed mesh
        viewer.data().clear();
        // Code for computing the mesh (V,F) from grid_points and grid_values
        if ((grid_points.rows() == 0) || (grid_values.rows() == 0)) {
            cerr << "Not enough data for Marching Cubes !" << endl;
            return true;
        }
        // Run marching cubes
        igl::copyleft::marching_cubes(grid_values, grid_points, resolution, resolution, resolution, V, F);
        if (V.rows() == 0) {
            cerr << "Marching Cubes failed!" << endl;
            return true;
        }

        igl::per_face_normals(V, F, FN);
        viewer.data().set_mesh(V, F);
        viewer.data().show_lines = true;
        viewer.data().show_faces = true;
        viewer.data().set_normals(FN);
    }

    return true;
}

bool callback_load_mesh(Viewer& viewer,string filename)
{
  igl::readOFF(filename,P,F,N);
  callback_key_down(viewer,'1',0);
  return true;
}

//helper boolean state for the interface 
bool state1 = true;
bool state2, state3;
int main(int argc, char *argv[]) {
    if (argc != 2) {
      cout << "Usage ex2_bin <mesh.off>" << endl;
      igl::readOFF("../data/bunny-1000.off",P,F,N);
    }
	else
	{
		// Read points and normals
		igl::readOFF(argv[1],P,F,N);
	}
    Viewer viewer;
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

	adjustAxis(); // alligning the axis

    viewer.callback_key_down = callback_key_down;

    menu.callback_draw_viewer_menu = [&]()
    {
      // Draw parent menu content
      menu.draw_viewer_menu();

      // Add new group
      if (ImGui::CollapsingHeader("Reconstruction Options", ImGuiTreeNodeFlags_DefaultOpen))
      {
        // Expose variable directly ...
        ImGui::InputInt("Resolution", &resolution, 0, 0);
        if (ImGui::Button("Reset Grid", ImVec2(-1,0)))
        {
          std::cout << "ResetGrid\n";
          // Recreate the grid
          createGrid();
          // Switch view to show the grid
          callback_key_down(viewer,'3',0);
        }

        // TODO: Add more parameters to tweak here...
		ImGui::SliderFloat("wendlandRadius", (float*)&wendlandRadius, 0.0, 5.0);

		ImGui::Text("Polynomial Degree:");
		if (ImGui::Checkbox("0 ", &state1)){
			polyDegree = 0;
			state2 = 0;
			state3 = 0;
		}
		ImGui::SameLine();
		if (ImGui::Checkbox("1 ", &state2))
		{
			polyDegree = 1;
			state1 = 0;
			state3 = 0;
		}
		ImGui::SameLine();
		if (ImGui::Checkbox("2 ", &state3))
		{
			polyDegree = 2;
			state1 = 0;
			state2 = 0;
		}
      }

    };

    callback_key_down(viewer, '1', 0);

    viewer.launch();
}
