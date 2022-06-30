#ifndef __DIJKSTRA__
#define __DIJKSTRA__

#include <igl/edges.h>
#include <igl/adjacency_list.h>
#include <cfloat>
#include "fibonacci.h"

void dijkstra(int p, std::vector<std::vector<int>> &adj_list, Eigen::MatrixXd &V, Eigen::VectorXd &Q, Eigen::VectorXi &prev);

void get_path(Eigen::MatrixXi prev, std::vector<std::vector<int>> &edges, int d);

#endif
