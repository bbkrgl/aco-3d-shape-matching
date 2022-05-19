#include "isocurve.h"

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

double f(double l)
{ return pow(l, 4); }

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

void distance_isocurve(int k, int i, int j, shape &I, shape &J)
{
	
}
