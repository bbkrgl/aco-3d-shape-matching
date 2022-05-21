#include "utils.h"
#include "dijkstra.h"

shape::shape(char* mesh_path)
{
	igl::readOFF(mesh_path, V, F);
	std::vector<std::vector<int>> adj_list;
	igl::adjacency_list(F, adj_list);

	dist_matrix = Eigen::MatrixXd(V.rows(), V.rows());
	for (int i = 0; i < V.rows(); i++) {
		Eigen::VectorXi prev(V.rows());
		Eigen::VectorXd Q(V.rows());
		dijkstra(i, adj_list, V, Q, prev);

		dist_matrix.row(i) = Q;
	}
}