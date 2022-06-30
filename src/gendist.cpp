#include "gendist.h"

void gendist_thread(int k, int start, int end, shape &I, shape &J, Eigen::MatrixXd &dist_matrix)
{
	for (int i = start; i < end; i++) {
		for (int j = 0; j < J.V.rows(); j++) {
			double d = distance_isocurve(k, i, j, I, J);
			std::cout << "Matrix " << "(" << i << ", " << j << ") = " << d << std::endl;
			dist_matrix(i, j) = d;
		}
	}
}

void gendist(int num_threads, int k, shape &I, shape &J, Eigen::MatrixXd &dist_matrix)
{
	int v_per_thread = I.V.rows() / num_threads;
	int start_i = 0;
	int end_i = v_per_thread;
	for (int i = 0; i < num_threads; i++) {
		//std::thread t(gendist_thread, k, start_i, end_i, I, J, dist_matrix);
		start_i = end_i;
		end_i += v_per_thread;
	}
}
