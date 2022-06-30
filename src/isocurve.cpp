#include "isocurve.h"
#define f(l) std::pow(l, 2)

int is_in_radius(double r, Eigen::MatrixXd &V, Eigen::VectorXi &face, Eigen::VectorXd &Q)
{ return (Q(face(0)) < r) + (Q(face(1)) < r) * 2 + (Q(face(2)) < r) * 4; }

double l_segment_len(double r, Eigen::MatrixXd &V, Eigen::VectorXi &face, Eigen::VectorXd &Q)
{
	int is_in_r = is_in_radius(r, V, face, Q);

	switch (is_in_r) {
	case 0b001:
	case 0b110:
	{
		double a1 = std::fabs((r - Q(face(0))) / (Q(face(1)) - Q(face(0))));
		double a2 = std::fabs((r - Q(face(0))) / (Q(face(2)) - Q(face(0))));

		Eigen::VectorXd v1 = (1 - a1) * V.row(face(0)) + a1 * V.row(face(1));
		Eigen::VectorXd v2 = (1 - a2) * V.row(face(0)) + a2 * V.row(face(2));
		
		return (v1 - v2).norm();
	}
	case 0b010:
	case 0b101:
	{
		double a1 = std::fabs((r - Q(face(1))) / (Q(face(0)) - Q(face(1))));
		double a2 = std::fabs((r - Q(face(1))) / (Q(face(2)) - Q(face(1))));

		Eigen::VectorXd v1 = (1 - a1) * V.row(face(1)) + a1 * V.row(face(0));
		Eigen::VectorXd v2 = (1 - a2) * V.row(face(1)) + a2 * V.row(face(2));

		return (v1 - v2).norm();
	}
	case 0b100:
	case 0b011:
	{
		double a1 = std::fabs((r - Q(face(2))) / (Q(face(0)) - Q(face(2))));
		double a2 = std::fabs((r - Q(face(2))) / (Q(face(1)) - Q(face(2))));

		Eigen::VectorXd v1 = (1 - a1) * V.row(face(2)) + a1 * V.row(face(0));
		Eigen::VectorXd v2 = (1 - a2) * V.row(face(2)) + a2 * V.row(face(1));

		return (v1 - v2).norm();
	}
	}
	return 0;
}

double isocurve_len_r(double r, Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXd &Q)
{
	double l = 0;
	for (int i = 0; i < F.rows(); i++) {
		Eigen::VectorXi face = F.row(i);
		l += l_segment_len(r, V, face, Q);
	}
	return l;
}

void isocurve(int p, int k, shape &S, Eigen::VectorXd &H)
{
	Eigen::VectorXd Q = S.dist_matrix.row(p);
	double max_len = Q.maxCoeff();

	double sc = f(max_len) / max_len;
	double inc = max_len / k;

	for (int i = 0; i < k; i++) {
		double r = f((i + 1) * inc) / sc;
		H(i) = isocurve_len_r(r, S.V, S.F, Q);
	}
}

double partial_distance_isocurve(int k, int i, int j, shape &I, shape &J)
{
	double dist = 0;

	// To avoid overflow, scale the distance matrix that will just be used to determine if a point is in the radius
	double max_len_I = I.dist_matrix.row(i).maxCoeff();
	double max_len_J = J.dist_matrix.row(j).maxCoeff();

	double min_dist = std::min(max_len_I, max_len_J);
	max_len_I = max_len_I / min_dist;
	max_len_J = max_len_J / min_dist;

	Eigen::VectorXd Q_I = I.dist_matrix.row(i) / min_dist;
	double sc_I = f(max_len_I) / max_len_I;

	Eigen::VectorXd Q_J = J.dist_matrix.row(j) / min_dist;

	for (int i = 0; i < k; i++) {
		double r_i = f((i + 1) * max_len_I / k) / sc_I; // TODO: Check r_i

		if (r_i > std::min(max_len_I, max_len_J))
			break;

		double d_i = isocurve_len_r(r_i, I.V, I.F, Q_I);
		double d_j = isocurve_len_r(r_i, J.V, J.F, Q_J);
		dist += std::abs((d_i - d_j) / (d_i + d_j));
	}

	return dist;
}

void initialize_dist_hists(int k, shape &I, shape &J)
{
	I.dist_hist = std::vector<std::vector<double>>(I.V.rows());
	J.dist_hist = std::vector<std::vector<double>>(J.V.rows());
	std::vector<double> r(k);

	std::cout << "Distance hist I" << std::endl;
	for (int i = 0; i < I.V.rows(); i++) {
		std::cout << "Point i " << i << std::endl;
		Eigen::VectorXd Q_I = I.dist_matrix.row(i);
		double max_len_I = Q_I.maxCoeff();
		double sc_I = f(max_len_I) / max_len_I;

		for (int i_ = 0; i_ < k; i_++) {
			double r_i = f((i_ + 1) * max_len_I / k) / sc_I; // TODO: Check r_i
			std::cout << "r_i " << r_i << std::endl;
			r[i_] = r_i;

			double d_i = isocurve_len_r(r_i, I.V, I.F, Q_I);
			std::cout << "d_i " << d_i << std::endl;
			I.dist_hist[i].push_back(d_i);
		}

		std::cout << "--------------------------" << std::endl;
		std::cout << std::endl;
	}

	std::cout << "Distance hist J" << std::endl;
	for (int j = 0; j < J.V.rows(); j++) {
		std::cout << "Point j " << j << std::endl;
		Eigen::VectorXd Q_J = J.dist_matrix.row(j);
		double max_len_J = Q_J.maxCoeff();

		for (int i_ = 0; i_ < k; i_++) {
			double r_i = r[i_];
			std::cout << "r_i " << r_i << std::endl;

			double d_j = isocurve_len_r(r_i, J.V, J.F, Q_J);
			std::cout << "d_j " << d_j << std::endl;
			J.dist_hist[j].push_back(d_j);
		}
		
		std::cout << "--------------------------" << std::endl;
		std::cout << std::endl;
	}
}

double distance_isocurve(int k, int i, int j, shape &I, shape &J)
{
	double d_ij = partial_distance_isocurve(k, i, j, I, J);
	double d_ji = partial_distance_isocurve(k, j, i, J, I);

	return (d_ij + d_ji) / 2;
}

double distance_isocurve_smart(int k, int i, int j, shape &I, shape &J)
{
	double d = 0;
	if (I.dist_hist.empty() || J.dist_hist.empty())
		initialize_dist_hists(k, I, J);

	for (int i_ = 0; i_ < k; i_++) {
		double r = std::abs(I.dist_hist[i][i_] - J.dist_hist[j][i_]) / (I.dist_hist[i][i_] + J.dist_hist[j][i_]);
		if (r != r)
			continue;
		d += std::abs(I.dist_hist[i][i_] - J.dist_hist[j][i_]) / (I.dist_hist[i][i_] + J.dist_hist[j][i_]);
	}

	return d;
}
