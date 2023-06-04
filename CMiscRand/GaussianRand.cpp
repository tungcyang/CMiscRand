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

static __m256i avxRand = _mm256_set_epi32(2000000000, 1000000000, 30000000, 1000000, 30000, 1000, 30, 0);
                                                                                // We use 0, 30, 1000 and 30000 as random seeds
                                                                                //     if the user does not call sGaussianRandVec().
static int	numAvailableResults = 0;

// sGaussianRandVec() resets the four random seeds used to generate Gaussian random variables for
//     GaussianRandVec().
void		__cdecl sGaussianRandVec(int s0, int s1, int s2, int s3,
	                                 int s4, int s5, int s6, int s7)
{
	avxRand = _mm256_set_epi32(s7, s6, s5, s4, s3, s2, s1, s0);
	numAvailableResults = 0;
}

// GaussianRandVec() returns a Gaussian-distributed double random number with zero mean and unit variance.
double		__cdecl	GaussianRandVec()
{
#define NUM_GAUSSIANRAND_GENERATED      16
	// The computations are based on AVX512, so we want results[] to be aligned on 64-byte.
	static __declspec(align(64)) double results[NUM_GAUSSIANRAND_GENERATED];

	if (numAvailableResults == 0)
	{
		// We don't have available results in results[]; time to refresh them in parallel.
		//     The implementation here assumes AVX512 support.  We carry out the computations in
		//     eight (NUM_GAUSSIANRAND_GENERATED/2) lanes of doubles (64-bit).  We only use 32-bit integers.
		const __m256i l214013 = _mm256_set_epi32(214013L, 214013L, 214013L, 214013L, 214013L, 214013L, 214013L, 214013L);
		const __m256i l2531011 = _mm256_set_epi32(2531011L, 2531011L, 2531011L, 2531011L, 2531011L, 2531011L, 2531011L, 2531011L);
		const __m256i m32767 = _mm256_set_epi32(0x7fff, 0x7fff, 0x7fff, 0x7fff, 0x7fff, 0x7fff, 0x7fff, 0x7fff);
		const __m512d dTwoOverRANDMAX = _mm512_set1_pd(2.0 / RAND_MAX);
		const __m512d dOnes = _mm512_set1_pd(1.0);
		const __m512d dZeros = _mm512_setzero_pd();
		const __m512d dMTwos = _mm512_set1_pd(-2.0);
		__m256i prevAvxRand = avxRand;
		__m512d avxdS = dZeros;
		__mmask8 dMasks = 0xff;
		__m512d avxdU, avxdV, dTmp;

		do
		{
			// dMasks is the mask for eight 64-bit lanes, while iMasks is the corresponding eight 32-bit lanes.
			// We use iMasks to choose if we want to update the avxRand lanes.  If the lane is 0x00000000, we
			//     restore avxRand lane to its prior value, and thus going through the loop just repeats.
			avxRand = _mm256_mask_blend_epi32(dMasks, prevAvxRand, avxRand);
			prevAvxRand = avxRand;

			// Evaluate avxRand = avxRand * 214013L + 2531011L and then (avxRand >> 16) & 0x7fff for u.
			//     Note that avxRand is between 0 and 32767, and avxRand * 214013L + 2531011L can
			//     overflow 32-bit integers.  However, since rand() internally only used int to keep its
			//     states, this means we can ignore the carry-outs.
			avxRand = _mm256_mullo_epi32(avxRand, l214013);				// avxRand in 8 32-bit integer lanes
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
			// Both avxdU and avxdV are in eight 64-bit double lanes.  Note that we use dMasks to decide if we want
			//     to update u and v.  Once dMasks shows 0 bit in one lane, the corresponding u and v are no longer
			//     updated, and this implies s also stops being updated.
			avxdU = _mm512_cvtepi32_pd(avxU);
			avxdV = _mm512_cvtepi32_pd(avxV);
			avxdU = _mm512_fmsub_pd(avxdU, dTwoOverRANDMAX, dOnes);
			avxdV = _mm512_fmsub_pd(avxdV, dTwoOverRANDMAX, dOnes);

			// s = u*u + v*v;
			avxdS = _mm512_add_pd(_mm512_mul_pd(avxdU, avxdU), _mm512_mul_pd(avxdV, avxdV));

			// Check if (s >= 1.0) || (s == 0.0).  Note that we evaluate the loop-terminating condition in
			//     eight lanes.  *ANY* lane with the specified condition means we continue the loop.  dMasks
			//     started with 0xff.  When one of the lanes hit 0 (represented by a bit in dMasks), it no
			//     longer needs to keep looping.  If dMasks becomes 0, the loop is broken.  Otherwise the
			//     loop continues, but the 0-bit lane will stop avxRand from being updated.
			dMasks = _mm512_cmp_pd_mask(avxdS, dOnes, _CMP_GE_OQ) | _mm512_cmp_pd_mask(avxdS, dZeros, _CMP_EQ_OQ);
		} while (dMasks);

		// u*sqrt(-2.0 * log(s)/s)
		// v*sqrt(-2.0 * log(s)/s)
		// It seems Intel Short Vector Math Library (SVML) has been supported by Microsoft Visual Studio
		//     2019 since Preview 2, thus we can use _mm256_log_pd() and _mm256_div_pd() here.
		//     https://devblogs.microsoft.com/cppblog/msvc-backend-updates-in-visual-studio-2019-preview-2/
		dTmp = _mm512_log_pd(avxdS);	// SVML/AVX
		dTmp = _mm512_div_pd(dTmp, avxdS);
		dTmp = _mm512_sqrt_pd(_mm512_mul_pd(dTmp, dMTwos));
		avxdU = _mm512_mul_pd(avxdU, dTmp);
		avxdV = _mm512_mul_pd(avxdV, dTmp);

		_mm512_store_pd(results, avxdU);
		_mm512_store_pd(results + NUM_GAUSSIANRAND_GENERATED/2, avxdV);

		numAvailableResults = NUM_GAUSSIANRAND_GENERATED;
	}

	numAvailableResults--;
	return results[numAvailableResults];
}
