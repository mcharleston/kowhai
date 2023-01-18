/*
 * approx.h
 *
 *  Created on: 18 Jan 2023
 *      Author: mac
 */

#ifndef UTILITY_APPROX_H_
#define UTILITY_APPROX_H_

#include <cmath>

class approx {
private:
	double tolerance;
public:
	approx(double t=1e-10) : tolerance(t) {}
	inline bool operator()(double d1, double d2) const {
		return (std::abs(d1-d2) <= tolerance);	// possible this could be sped up cleverly but I don't want to yet.
	}
};


#endif /* UTILITY_APPROX_H_ */
