#ifndef __ISOCURVE__
#define __ISOCURVE__

#include "utils.h"

void isocurve(int p, int k, shape &S, Eigen::VectorXd &H);
double distance_isocurve(int k, int i, int j, shape &I, shape &J);
void initialize_dist_hists(int k, shape &I, shape &J);
double distance_isocurve_smart(int k, int i, int j, shape &I, shape &J);

#endif
