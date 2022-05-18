#include "aco.h"
#include <igl/opengl/glfw/Viewer.h>


int main(int argc, char *argv[])
{
	if (argc < 3) {
		std::cerr << "Usage: ./3daco display <mesh_path>" << std::endl;
		std::cerr << "Usage: ./3daco <mesh1_path> <mesh2_path>" << std::endl;
		return 0;
	}

	if (!strcmp("display", argv[1])) {
		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		igl::readOFF(argv[2], V, F);

		igl::opengl::glfw::Viewer viewer;
		viewer.data().set_mesh(V, F);
		viewer.data().set_face_based(true);
		viewer.launch();
	} else {
		shape I(argv[1]);
		shape J(argv[2]);
		J.V.col(0).array() += 100;

		best_mc m;
		//shape_match(I, J, m);

		Eigen::MatrixXd V(I.V.rows() + J.V.rows(), I.V.cols());
		V << I.V, J.V;
		Eigen::MatrixXi F(I.F.rows() + J.F.rows(), I.F.cols());
		F << I.F, (J.F.array() + I.V.rows());

		Eigen::MatrixXd C(F.rows(), 3);
		C << Eigen::RowVector3d(0.2, 0.3, 0.8).replicate(I.F.rows(), 1),
			Eigen::RowVector3d(1.0, 0.7, 0.2).replicate(J.F.rows(), 1);

		igl::opengl::glfw::Viewer viewer;
		viewer.data().set_mesh(V, F);
		viewer.data().set_colors(C);
		viewer.data().set_face_based(true);
		viewer.launch();
	}
}