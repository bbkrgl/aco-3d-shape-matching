#ifndef __ACO3D__
#define __ACO3D__

#include <pagmo/algorithm.hpp>
#include <pagmo/algorithms/gaco.hpp>

struct aco3d_problem {
	pagmo::vector_double fitness(const pagmo::vector_double &dv) const;
	std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const;
};

#endif
