#include "utils.h"
#include <string>
#include <fstream>

shape::shape(char* mesh_path)
{
	std::cout << "Geodesic distances of " << mesh_path << std::endl;
	igl::readOFF(mesh_path, V, F);
	igl::adjacency_list(F, this->adj_list);

	dist_matrix = Eigen::MatrixXd(V.rows(), V.rows());
	for (int i = 0; i < V.rows(); i++) {
		Eigen::VectorXi prev(adj_list.size());
		Eigen::VectorXd Q(adj_list.size());
		dijkstra(i, this->adj_list, V, Q, prev);
		dist_matrix.row(i) = Q;
	}

	std::string filename = (char*) (std::strrchr(mesh_path, '/') + 1);
	std::ofstream file(filename + "_dist");
	if (file.is_open()) {
		file << dist_matrix;
		file.close();
	}
}
