#ifndef __ACO3D__
#define __ACO3D__

#include "utils.h"

typedef struct graph {

	std::vector<std::vector<double>> pheromones;

	graph() {};
	graph(shape &I, shape &J);
	
	matching construct_matching();
	void update_pheromones();
} graph;

void shape_match(shape &I, shape &J, best_mc &m);

#endif
