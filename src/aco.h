#ifndef __ACO3D__
#define __ACO3D__

#include <vector>
#include <igl/readOFF.h>
#include <dlib/optimization/max_cost_assignment.h>

typedef std::pair<int, int> matching;
typedef std::pair<matching, double> best_mc;

typedef struct shape {
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	shape() {};
	shape(char* mesh_path) { igl::readOFF(mesh_path, V, F); };
} shape;

typedef struct graph {

	std::vector<std::vector<double>> pheromones;

	graph() {};
	graph(shape &I, shape &J);
	
	matching construct_matching();
	void update_pheromones();
} graph;

void shape_match(shape &I, shape &J, best_mc &m);

#endif