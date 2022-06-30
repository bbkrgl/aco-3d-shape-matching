#include <igl/opengl/glfw/Viewer.h>

#include "aco.h"
#include "isocurve.h"

void write_matching_to_file(std::vector<int>& m)
{
	std::ofstream fout("matching.txt");
	for (int i = 0; i < m.size(); i++)
		fout << m[i] << std::endl;
	fout.close();
}

void read_matchings_from_file(std::vector<int>& m)
{
	std::ifstream fin("matching.txt");
	while (!fin.eof()) {
		int m_;
		fin >> m_;
		m.push_back(m_);
	}
	fin.close();
}

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
		
		Eigen::MatrixXd C_I(F.rows(), 3);
		C_I << Eigen::RowVector3d(0.2, 0.3, 0.8).replicate(F.rows(), 1);

		igl::opengl::glfw::Viewer viewer;
		viewer.data().set_mesh(V, F);
		viewer.data().set_colors(C_I);
		viewer.data().set_face_based(true);
		viewer.launch();
	} else {
		shape I(argv[1]);
		shape J(argv[2]);
		J.V.col(0).array() += J.V.col(0).maxCoeff()*10;

		int n = atoi(argv[3]);
		//std::vector<int> matching;
		//hungarian_matching(n, I, J, matching);

		aco a(I, J, argv[4]);
		a.match();

		std::vector<int>* matching = a.best_matching;
		//write_matching_to_file(*matching);

		Eigen::MatrixXd C_I(I.F.rows(), 3);
		C_I << Eigen::RowVector3d(0.2, 0.3, 0.8).replicate(I.F.rows(), 1); // Blue

		Eigen::MatrixXd C_J(J.F.rows(), 3);
		C_J << Eigen::RowVector3d(1.0, 0.7, 0.2).replicate(J.F.rows(), 1); // Yellow

		igl::opengl::glfw::Viewer viewer;
		viewer.data().set_mesh(I.V, I.F);
		viewer.data().set_colors(C_I);

		int nm = viewer.append_mesh();
		viewer.data(nm).set_mesh(J.V, J.F);
		viewer.data(nm).set_colors(C_J);
	
		Eigen::MatrixXd matching_I(n, 3);
		Eigen::MatrixXd matching_J(n, 3);		
		for (int i = 0; i < n; i++) {
			matching_I.row(i) = I.V.row(i);
			matching_J.row(i) = J.V.row((*matching)[i]);
		}
		viewer.data().add_edges(matching_I, matching_J, Eigen::RowVector3d(1, 0, 0));

		viewer.launch();
	}
}
