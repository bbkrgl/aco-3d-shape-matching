#include "utils.h"
#include "dijkstra.h"

shape::shape(char* mesh_path)
{
	igl::readOFF(mesh_path, V, F);
	igl::adjacency_list(F, this->adj_list);

	dist_matrix = Eigen::MatrixXd(V.rows(), V.rows());
	for (int i = 0; i < V.rows(); i++) {
		Eigen::VectorXi prev(V.rows());
		Eigen::VectorXd Q(V.rows());
		dijkstra(i, this->adj_list, V, Q, prev);

		dist_matrix.row(i) = Q;
	}
}