#ifndef __ACO3D__
#define __ACO3D__

#include <cmath>
#include <iostream>
#include <random>
#include "Hungarian.h"
#include "utils.h"
#include "isocurve.h"
#include <fstream>

#define NUM_ITERS 500
#define NUM_ANTS 1
#define k 8
#define v 0.7
#define alpha 0.3
#define phi 0.1
#define delta 0.01
#define th0 1

typedef struct graph {
	shape &I;
	shape &J;

	double th_min;
	double s_R;
	double s_I;

	Eigen::MatrixXd ph_I_J;
	Eigen::MatrixXd dist_matrix;

	Eigen::MatrixXd dist_matrix_I;
	Eigen::MatrixXd dist_matrix_J;

	graph(shape &I, shape &J): I(I), J(J)
	{
		dist_matrix_I = I.dist_matrix.rowwise().normalized();
		dist_matrix_J = J.dist_matrix.rowwise().normalized();

		ph_I_J = Eigen::MatrixXd(I.V.rows(), J.V.rows());
		ph_I_J.fill(th0);
		dist_matrix = Eigen::MatrixXd(I.V.rows(), J.V.rows());
	}

	graph(shape &I, shape &J, char* dist_matrix_name): I(I), J(J)
	{	
		dist_matrix_I = I.dist_matrix.rowwise().normalized();
		dist_matrix_J = J.dist_matrix.rowwise().normalized();

		ph_I_J = Eigen::MatrixXd(I.V.rows(), J.V.rows());
		ph_I_J.fill(th0);
		dist_matrix = Eigen::MatrixXd(I.V.rows(), J.V.rows());

		std::ifstream fin(dist_matrix_name, std::ios::in);
		if (/*fin.is_open()*/0) { // TODO: Fix file read
			for (int i = 0; i < I.V.rows(); i++) {
				for (int j = 0; j < J.V.rows(); j++) {
					double d = 0;
					fin >> d;
					dist_matrix(i, j) = d;
					std::cerr << strerror(errno) << std::endl;
				}
			}
			fin.close();
		} else {
			std::cout << "-----------------------------" << std::endl;
			std::cout << "Isocurves:" << std::endl;
			for (int i = 0; i < I.V.rows(); i++) {
				for (int j = 0; j < J.V.rows(); j++) {
					double d = distance_isocurve_smart(k, i, j, I, J);
					std::cout << "Matrix " << "(" << i << ", " << j << ") = " << d << std::endl;
					dist_matrix(i, j) = d;
				}
			}

			dist_matrix.rowwise().normalize();

			std::ofstream fout(dist_matrix_name);
			if (fout.is_open()) {
				fout << dist_matrix;
				fout.close();
			}
			std::cout << "-----------------------------" << std::endl;
		}
		
		th_min = 0.1 / I.V.rows();
		s_R = 0.1 * dist_matrix_I.maxCoeff();
		s_I = 0.1 * dist_matrix.maxCoeff();
	}

	void update_pheromones(std::vector<std::pair<std::vector<int>*, double>> matchings);
} graph;

typedef struct aco {
	graph G;
	std::vector<std::vector<int>*> matchings;
	std::vector<int>* best_matching = 0;
	double best_matching_cost = DBL_MAX;

	aco(shape &I, shape &J, char* dist_matrix_name): G(I, J, dist_matrix_name) {}

	~aco()
	{
		for (int i = 0; i < matchings.size(); i++)
			delete matchings[i];
	}

	double edge_probability(std::vector<int> &m, int i, int j, int i_, int i__);
	double n_ij(std::vector<int> &m, int i, int j, int i_, int i__);
	std::vector<int>* construct_matching();
	void match();
} aco;

double similarity_term(graph &G, std::vector<int> &m);
double proximity_term(graph &G, std::vector<int>& m);
double qap(graph &G, std::vector<int>& m);

void hungarian_matching(int max_i, shape &I, shape &J, std::vector<int> &m);

#endif
