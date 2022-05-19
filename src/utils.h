#ifndef __UTIL_ACO3D__
#define __UTIL_ACO3D__

#include <vector>
#include <igl/readOFF.h>
#include "dijkstra.h"
#include "igl/adjacency_list.h"

typedef std::pair<int, int> matching;
typedef std::pair<matching, double> best_mc;

typedef struct shape {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::MatrixXd dist_matrix;

	shape() {};
	shape(char* mesh_path);
} shape;

#endif
