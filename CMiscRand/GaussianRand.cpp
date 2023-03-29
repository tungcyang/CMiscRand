// GaussianRand.cpp implements the Polar form of Box–Muller transform (https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform)
// to generate a pair of independent zero-mean, unit-variance Gaussian-distributed random variables based on a pair of
// independent uniform distribution variables.

#include <stdlib.h>
#include <math.h>
#include <immintrin.h>
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

static __m256i avxRand = _mm256_set_epi32(0, 0, 0, 0, 30000, 1000, 30, 0);      // We use 0, 30, 1000 and 30000 as random seeds
                                                                                //     if the user does not call sGaussianRandVec().
static int	numAvailableResults = 0;

// sGaussianRandVec() resets the four random seeds used to generate Gaussian random variables for
//     GaussianRandVec().
void		__cdecl sGaussianRandVec(int s0, int s1, int s2, int s3)
{
	avxRand = _mm256_set_epi32(0, 0, 0, 0, s3, s2, s1, s0);
	numAvailableResults = 0;
}

// GaussianRandVec() returns a Gaussian-distributed double random number with zero mean and unit variance.
double		__cdecl	GaussianRandVec()
{
#define NUM_GAUSSIANRAND_GENERATED  8
	// The computations are based on AVX/AVX2, so we want results[] to be aligned on 32-byte.
	static __declspec(align(32)) double results[NUM_GAUSSIANRAND_GENERATED];

	if (numAvailableResults == 0)
	{
		// We don't have available results in results[]; time to refresh them in parallel.
		//     The implementation here assumes AVX/AVX2 support.  We carry out the computations in
		//     four (NUM_GAUSSIANRAND_GENERATED/2) lanes of doubles (64-bit).  We only use 32-bit integers.
		const __m256i l214013 = _mm256_set_epi32(0, 0, 0, 0, 214013L, 214013L, 214013L, 214013L);
		const __m256i l2531011 = _mm256_set_epi32(0, 0, 0, 0, 2531011L, 2531011L, 2531011L, 2531011L);
		const __m256i m32767 = _mm256_set_epi32(0, 0, 0, 0, 0x7fff, 0x7fff, 0x7fff, 0x7fff);
		const __m256i bPermuteSelect = _mm256_set_epi32(0, 0, 0, 0, 6, 4, 2, 0);
		const __m256d dTwoOverRANDMAX = _mm256_set1_pd(2.0 / RAND_MAX);
		const __m256d dOnes = _mm256_set1_pd(1.0);
		const __m256d dZeros = _mm256_setzero_pd();
		const __m256d dMTwos = _mm256_set1_pd(-2.0);
		__m256i prevAvxRand = avxRand;
		__m256d avxdS = dZeros;
		__m256d dMasks = _mm256_castsi256_pd(_mm256_set1_epi64x(0xffffffffffffffff));
		__m256d avxdU, avxdV, dTmp;

		do
		{
			__m256i iMasks = _mm256_permutevar8x32_epi32(_mm256_castpd_si256(dMasks), bPermuteSelect);
			// dMasks is the mask for four 64-bit lanes, while iMasks is the corresponding four 32-bit lanes in the lower 128 bits.
			// We use iMasks to choose if we want to update the avxRand lanes.  If the lane is 0x00000000, we
			//     restore avxRand lane to its prior value, and thus going through the loop just repeats.
			avxRand = _mm256_or_si256(_mm256_and_si256(iMasks, avxRand),
				                      _mm256_andnot_si256(iMasks, prevAvxRand));
			prevAvxRand = avxRand;

			// Evaluate avxRand = avxRand * 214013L + 2531011L and then (avxRand >> 16) & 0x7fff for u.
			//     Note that avxRand is between 0 and 32767, and avxRand * 214013L + 2531011L can
			//     overflow 32-bit integers.  However, since rand() internally only used int to keep its
			//     states, this means we can ignore the carry-outs.
			avxRand = _mm256_mullo_epi32(avxRand, l214013);    // avxRand in 4 32-bit integer lanes
			avxRand = _mm256_add_epi32(avxRand, l2531011);
			__m256i     avxU = _mm256_srli_epi32(avxRand, 16);
			avxU = _mm256_and_si256(avxU, m32767);

			// Repeating the steps to generate v
			avxRand = _mm256_mullo_epi32(avxRand, l214013);
			avxRand = _mm256_add_epi32(avxRand, l2531011);
			__m256i     avxV = _mm256_srli_epi32(avxRand, 16);
			avxV = _mm256_and_si256(avxV, m32767);

			// u = 2.0*((double)rand() / RAND_MAX) - 1;
			// v = 2.0*((double)rand() / RAND_MAX) - 1;
			// Both avxdU and avxdV are in four 64-bit double lanes.  Note that we use dMasks to decide if we want
			//     to update u and v.  Once dMasks shows 0x0000 in one lane, the corresponding u and v are no longer
			//     updated, and this implies s also stops being updated.
			avxdU = _mm256_cvtepi32_pd(_mm256_castsi256_si128(avxU));
			avxdV = _mm256_cvtepi32_pd(_mm256_castsi256_si128(avxV));
			avxdU = _mm256_mul_pd(avxdU, dTwoOverRANDMAX);
			avxdV = _mm256_mul_pd(avxdV, dTwoOverRANDMAX);
			avxdU = _mm256_sub_pd(avxdU, dOnes);		// Using _mm256_fmadd_pd() here requires CPU FMA flag
			avxdV = _mm256_sub_pd(avxdV, dOnes);

			// s = u*u + v*v;
			avxdS = _mm256_add_pd(_mm256_mul_pd(avxdU, avxdU), _mm256_mul_pd(avxdV, avxdV));

			// Check if (s >= 1.0) || (s == 0.0).  Note that we evaluate the loop-terminating condition in
			//     four lanes.  *ANY* lane with the specified condition means we continue the loop.  dMasks started
			//     with four lanes of 0xffffffff.  When one of the lanes hit 0x00000000, it no longer needs to keep
			//     looping.  If all lanes give 0x00000000, _mm256_testz_pd() will return 1 and the loop is broken.
			//     Otherwise the loop continues, but the 0x00000000 lane will stop avxRand from being updated.
			dMasks = _mm256_or_pd(_mm256_cmp_pd(avxdS, dOnes, _CMP_GE_OQ), _mm256_cmp_pd(avxdS, dZeros, _CMP_EQ_OQ));
				// ordered non-signaling comparison
		} while (!_mm256_testz_pd(dMasks, dMasks));

		// u*sqrt(-2.0 * log(s)/s)
		// v*sqrt(-2.0 * log(s)/s)
		// Unfortunately there is no intrinsic for logarithms nor divisions, and there is no intrinsic for
		//     reciprocals in double.
		_mm256_store_pd(results, avxdS);
		for (int i = 0; i < NUM_GAUSSIANRAND_GENERATED/2; i++)
			results[i] = log(results[i])/results[i];
		dTmp = _mm256_load_pd(results);     // log(s)/s
		dTmp = _mm256_sqrt_pd(_mm256_mul_pd(dTmp, dMTwos));
		avxdU = _mm256_mul_pd(avxdU, dTmp);
		avxdV = _mm256_mul_pd(avxdV, dTmp);

		_mm256_store_pd(results, avxdU);
		_mm256_store_pd(results + NUM_GAUSSIANRAND_GENERATED/2, avxdV);

		numAvailableResults = NUM_GAUSSIANRAND_GENERATED;
	}

	numAvailableResults--;
	return results[numAvailableResults];
}
