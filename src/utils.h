#ifndef __UTIL_ACO3D__
#define __UTIL_ACO3D__

#include <vector>
#include <igl/readOFF.h>
#include "dijkstra.h"
#include "igl/adjacency_list.h"

typedef struct shape {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	std::vector<std::vector<int>> adj_list;
	Eigen::MatrixXd dist_matrix;

	shape() {};
	shape(char* mesh_path);
} shape;

#endif