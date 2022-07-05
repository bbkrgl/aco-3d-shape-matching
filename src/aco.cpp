#include "aco.h"

std::random_device rd;
std::mt19937 gen(rd());

int get_rand_int(std::vector<int>& available)
{
	std::uniform_int_distribution<> dis(0, available.size() - 1);
	int i = dis(gen);
	int ret = available[i];

	available.erase(available.begin() + i);

	return ret;
}

void graph::update_pheromones(std::vector<std::pair<std::vector<int>*, double>> matchings)
{
	for (int j = 0; j < ph_I_J.cols(); j++) {
		for (int i = 0; i < ph_I_J.rows(); i++) {
			ph_I_J(i, j) *= (1 - phi);
			if (ph_I_J(i, j) < th_min)
				ph_I_J(i, j) = th_min;
		}
	}

	for (int i = 0; i < matchings.size(); i++) {
		std::vector<int>* m = matchings[i].first;
		for (int i_ = 0; i_ < m->size(); i_++)
			ph_I_J(i_, (*m)[i_]) += delta / matchings[i].second; // TODO: Check
	}
}

double aco::n_ij(std::vector<int> &m, int i, int j, int i_, int i__)
{
	return std::exp(-std::pow(G.dist_matrix(i, j), 2) / this->G.s_R) * 
		(1 - std::exp(-std::pow(this->G.dist_matrix_I(i, i_), 2) / this->G.s_I) * std::abs(this->G.dist_matrix_I(i, i_) - this->G.dist_matrix_J(j, m[i_]))) * 
		(1 - std::exp(-std::pow(this->G.dist_matrix_I(i, i__), 2) / this->G.s_I) * std::abs(this->G.dist_matrix_I(i, i__) - this->G.dist_matrix_J(j, m[i__])));
}

double aco::edge_probability(std::vector<int> &m, int i, int j, int i_, int i__)
{
	double denom = 0;
	for (int l = 0; l < this->G.J.V.rows(); l++) // TODO: Figure out if N_i == J
		denom += (alpha * this->G.ph_I_J(i, l) + (1 - alpha) * n_ij(m, i, l, i_, i__));
	return (alpha * this->G.ph_I_J(i, j) + (1 - alpha) * n_ij(m, i, j, i_, i__)) / denom;
}

std::vector<int>* aco::construct_matching()
{
	std::vector<int> v_available(this->G.I.V.rows());
	for (int i = 0; i < v_available.size(); i++)
		v_available[i] = i;

	std::vector<int>* matching = new std::vector<int>(this->G.I.V.rows(), 0);
	int i = get_rand_int(v_available);

	std::vector<double> P = std::vector<double>(G.J.V.rows(), 0);
	int i_ = 0;
	int i__ = 0;
	while (!v_available.empty()) {
		//std::cout << "Current vertex: " << i << std::endl;
		//std::cout << "Current progress: " << G.I.V.rows() - v_available.size() << "/" << G.I.V.rows() << std::endl;

		double sum = 0;
		for (int j = 0; j < G.J.V.rows(); j++) // TODO: Figure out order preservation
			P[j] = edge_probability(*matching, i, j, i_, i__);

		std::discrete_distribution<> dis(P.begin(), P.end());
		(*matching)[i] = dis(gen);
		i = get_rand_int(v_available);
		i__ = i_;
		i_ = i;
	}

	return matching;
}

void aco::match()
{
	std::cout << "-----------------------------" << std::endl;
	std::cout << "Start Matching" << std::endl;

	for (int i = 0; i < NUM_ITERS; i++) {
		std::vector<std::pair<std::vector<int>*, double>> matchings;
		for (int j = 0; j < NUM_ANTS; j++) {
			std::cout << "Iteration " << i << " Ant " << j << std::endl;

			std::vector<int> *matching = construct_matching();
			double cost = qap(this->G, *matching);
			std::cout << "Cost " << cost << std::endl;
			matchings.push_back(std::make_pair(matching, cost));

			if (this->best_matching_cost > cost) {
				this->best_matching_cost = cost;
				this->best_matching = matching;
			}
		}
		this->G.update_pheromones(matchings);
	}

	std::cout << "Matching end" << std::endl;
	std::cout << "-----------------------------" << std::endl;
}

void hungarian_matching(int max_i, shape &I, shape &J, std::vector<int> &m)
{
	HungarianAlgorithm alg;
	std::vector<std::vector<double>> cost;
	std::cout << "Start cost computation" << std::endl;

	for (int i = 0; i < max_i; i++) {
		std::vector<double> c;
		for (int j = 0; j < J.V.rows(); j++) {
			double c_ = distance_isocurve_smart(k, i, j, I, J);
			c.push_back(c_);
			std::cout << "Cost: " << c_ << std::endl;
		}
		cost.push_back(c);
	}

	std::cout << "Start assignment computation" << std::endl;
	alg.Solve(cost, m);
}

double similarity_term(graph &G, std::vector<int> &m)
{
	double res = 0;
	for (int i = 0; i < G.I.V.rows(); i++)
		res += std::abs(1 - std::exp(-std::pow(G.dist_matrix(i, m[i]), 2) / G.s_R));

	return res / G.I.V.rows();
}

double proximity_term(graph &G, std::vector<int> &m)
{
	double res = 0;
	for (int i = 0; i < G.I.V.rows(); i++) {
		for (int i_ = 0; i_ < G.I.V.rows(); i_++) {
			if (i == i_)
				continue;
			res += std::exp(-std::pow(G.dist_matrix_I(i, i_), 2) / G.s_I)
				* std::abs(G.dist_matrix_I(i, i_) - G.dist_matrix_J(m[i], m[i_]));
		}
	}

	return res / (G.I.V.rows() * (G.I.V.rows() - 1)) / 2;
}

double qap(graph &G, std::vector<int> &m)
{
	double sim = similarity_term(G, m);
	double prox = proximity_term(G, m);

	std::cout << "Similarity: " << sim << std::endl;
	std::cout << "Proximity: " << prox << std::endl;

	return (1 - v) * sim + v * prox;
}
