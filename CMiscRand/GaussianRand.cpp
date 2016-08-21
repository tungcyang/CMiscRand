// GaussianRand.cpp implements the Polar form of Box–Muller transform (https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform)
// to generate a pair of independent zero-mean, unit-variance Gaussian-distributed random variables based on a pair of
// independent uniform distribution variables.

#include <stdlib.h>
#include <math.h>
#include "MiscRand.h"

// GaussianRand() returns a Gaussian-distributed double random number with zero mean and unit variance.
double		__cdecl	GaussianRand()
{
	double	u, v, s;

	do {
		// Evaluating u and v, where they are independent and uniformly distributed in the closed interval [?1, +1], 
		u = 2.0 * ((double) rand()/RAND_MAX) - 1;
		v = 2.0 * ((double) rand()/RAND_MAX) - 1;

		// Evaluating s = R*R = u*u + v*v.  Discarding s if s == 0.0 or s >= 1.0.
		s = u*u + v*v;
	} while ((s >= 1.0) || (s == 0.0));

	// We are sure s is not 0.0 here.  Note that we evaluate z0 only here.
	return (u*sqrt(-2.0 * log(s)/s));
}