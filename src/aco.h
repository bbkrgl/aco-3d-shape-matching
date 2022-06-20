#ifndef __ACO3D__
#define __ACO3D__

#include <cmath>
#include <iostream>
#include <random>
#include "Hungarian.h"


#include "utils.h"
#include "isocurve.h"

#define k 30
#define v 3
#define alpha 0.3
#define phi 0.1
#define s_R 0.1
#define s_I 0.1
#define th0 0.1
#define delta 0.1
#define th_min 0.01
#define NUM_ITERS 1000
#define NUM_ANTS 1

typedef struct graph {
	shape &I;
	shape &J;

	Eigen::MatrixXd ph_I_J;

	graph(shape &I, shape &J): I(I), J(J)
	{
		ph_I_J = Eigen::MatrixXd(I.V.rows(), J.V.rows());
		ph_I_J.fill(th0);
	};

	void update_pheromones(std::vector<std::pair<std::vector<int>*, double>> matchings);
} graph;

typedef struct aco {
	graph G;
	std::vector<std::vector<int>*> matchings;
	std::vector<int>* best_matching = 0;
	double best_matching_cost = DBL_MAX;

	aco(shape &I, shape &J): G(I, J) {};

	~aco()
	{
		for (int i = 0; i < matchings.size(); i++)
			delete matchings[i];
	};

	double edge_probability(std::vector<int> &m, int i, int j, int i_, int i__);
	double n_ij(std::vector<int> &m, int i, int j, int i_, int i__);
	std::vector<int>* construct_matching();
	void match();
} aco;

double similarity_term(std::vector<int>& m, shape &I, shape &J);
double proximity_term(std::vector<int>& m, shape &I, shape &J);
double qap(std::vector<int>& m, shape &I, shape &J);

void hungarian_matching(int max_i, shape &I, shape &J, std::vector<int> &m);

#endif
