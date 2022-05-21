#ifndef __ISOCURVE__
#define __ISOCURVE__

#include "utils.h"

void isocurve(int p, int k, shape &S, Eigen::VectorXd &H);
double distance_isocurve(int i, int j, shape &I, shape &J);

#endif