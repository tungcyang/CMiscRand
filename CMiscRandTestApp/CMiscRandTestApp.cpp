// CMiscRandTestApp.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <stdlib.h>
#include <immintrin.h>
#include <intrin.h>
#include "MiscRand.h"
#include "cpudetect.h"

#define	RAND_SEED			0
#define RAND_SEQUENCE_LEN	8
#define NUM_GAUSSIAN_RANDOM		524288
#define GAUSSIAN_BIN_RANGE		4			// We try to put the generated Gaussian random numbers into bins in
											// [-GAUSSIAN_BIN_RANGE, GAUSSIAN_BIN_RANGE]
#define NUM_GAUSSIAN_BINS		32
#define NUM_RDRAND_ITERATIONS	1000000
#define RDRAND_RETRY            10

int main(int argc, char* argv[])
{
	unsigned long   i;
#if 0
	// In the first part, we evaluate how soon rand() seemingly repeats its returned value.  Note
	// that this does not necessarily equal the period of the pseudo-random sequence rand() generates.
	{
		unsigned long	uCount;
		int		iFirstRandVal, iRandVal;

		fprintf(stdout, "\n\n====== Part One ======\n");

		fprintf(stdout, "RAND_MAX is %d\n", RAND_MAX);
		srand(RAND_SEED);
		iFirstRandVal = rand();
		uCount = 0;
		while (true)
		{
			iRandVal = rand();
			uCount++;
			if (iRandVal == iFirstRandVal)
				break;
		}
		fprintf(stdout, "With %d as srand() seed, the first repeated rand() value is %d\nafter %d calls.", RAND_SEED, iRandVal, uCount);
	}

	// In the second part, we attempt to find the period of the pseudo-random sequence rand() generates
	// by comparing the partial sequences separated by the period.  The longer the sequences, the more
	// likely we found the period.  Note that part one is a special case for sequences of length 1.
	{
		unsigned long	uCount;
		int		iInitRandSequence[RAND_SEQUENCE_LEN], iRandValSequence[RAND_SEQUENCE_LEN];

		fprintf(stdout, "\n\n====== Part Two ======\n");

		// Filling the initial subsequence from rand().
		srand(RAND_SEED);
		for (i = 0; i < RAND_SEQUENCE_LEN; i++)
			iInitRandSequence[i] = rand();

		// Note that RAND_MAX is 32767, but rand() period is much larger, 2147483648,
		// by Microsoft Visual Studio Community 2015.
		// It is believed the following expression is used for computing rand():
		//     ((holdrand = holdrand * 214013L + 2531011L) >> 16) & 0x7fff;
		// "rand() questions", in http://www.programmersheaven.com/discussion/382052/rand-questions

		// Starting the loop to determine the period of pseudo-random number generator rand().
		uCount = RAND_SEQUENCE_LEN;	// We have already called rand() RAND_SEQUENCE_LEN times.
		while (1)
		{
			// Retrieve the next rand().  If it is different from iFirstRand[0], we do not have a
			// repeat pattern yet and we can jump to the next immediately.
			for (i = 0; i < RAND_SEQUENCE_LEN; i++)
			{
				iRandValSequence[i] = rand();
				uCount++;

				if (iRandValSequence[i] != iInitRandSequence[i])
					break;
			}

			// There are two possibilities for us to be here:
			// (1) i is RAND_SEQUENCE_LEN -- we have a repetition pattern RAND_SEQUENCE_LEN rand() calls before,
			//     so subtract RAND_SEQUENCE_LEN from uCount and we are done.
			// (2) i < RAND_SEQUENCE_LEN and iRandValSequence[i] != iInitRandSequence[i] -- we don't have a
			//     repetition pattern yet.  Keep going.
			if (i == RAND_SEQUENCE_LEN)
			{
				uCount -= RAND_SEQUENCE_LEN;
				break;
			}
		}
		fprintf(stdout, "rand() period is %u\n", uCount);

		// In the third part, we try to check if LargerRand() indeed generates the same pseudo-random
		// numbers as rand() before rand() repeats if (long) LargerRand() is mapped to int appropriately.
		int		iRandVal;
		int		iLargerRandVal;

		fprintf(stdout, "\n\n====== Part Three ======\n");

		srand(RAND_SEED);
		sLargerRand(RAND_SEED);
		for (i = 0; i < uCount; i++)
		{
			iRandVal = rand();
			iLargerRandVal = (int)((LargerRand() >> 16) & 0x7fff);

			if (iRandVal != iLargerRandVal)
			{
				fprintf(stdout, "rand() and ((LargerRand() >> 16) & 0x7fff) return different values: %d and %d\n",
					iRandVal, iLargerRandVal);
				break;
			}
		}
		if (i == uCount)
		{
			fprintf(stdout, "rand() and ((LargerRand() >> 16) & 0x7fff) return identical values for all %u calls.\n",
				uCount);
		}
	}

	// In the fourth part, we call the Gaussian random number generator many times and have a rough idea of the
	// "Bell curve".
	{
		unsigned int	uGaussianBinCounts[NUM_GAUSSIAN_BINS];
		double			dfGaussianBinGap = 2.0*GAUSSIAN_BIN_RANGE / NUM_GAUSSIAN_BINS;

		// Initialization
		fprintf(stdout, "\n\n====== Part Four ======\n");
		srand(0);
		for (i = 0; i < NUM_GAUSSIAN_BINS; i++)
			uGaussianBinCounts[i] = 0;

		// Generating the Gaussian random numbers.
		for (i = 0; i < NUM_GAUSSIAN_RANDOM; i++)
		{
			double	dfGaussianRandNum = GaussianRand();

			// Which bin will dfGaussianRandNum fall in?
			int		iBinIndex = (int) ((dfGaussianRandNum + GAUSSIAN_BIN_RANGE)/dfGaussianBinGap);
			if ((iBinIndex >= 0) && (iBinIndex < NUM_GAUSSIAN_BINS))
				uGaussianBinCounts[iBinIndex]++;
		}

		// Printing out the Gaussian bin counts.
		for (i = 0; i < NUM_GAUSSIAN_BINS; i++)
			fprintf(stdout, "%2dth bin (from %+f to %+f) contains %5d samples.\n", i,
				(-GAUSSIAN_BIN_RANGE + i * dfGaussianBinGap), (-GAUSSIAN_BIN_RANGE + (i + 1) * dfGaussianBinGap),
				uGaussianBinCounts[i]);
	}

	// In the fifth part, we generate many uniformly distributed random variables in [0, 1) with Intel/AMD
	//     RDRAND instruction.
	{
		unsigned int      uRdSample;
		double            fRdSample;
		double            fSumSamples = 0.0, fSumSampleSquares = 0.0;
		double            fSumCycles = 0.0, fNumRetries = 0.0;
		unsigned int      uRetryCount, uNumGoodSamples, uMaxRetries = 0;
		unsigned int      uMinRDCycles = UINT_MAX, uMinRDCount = 0;
		unsigned int      uMaxRDCycles = 0, uMaxRDCount = 0;
		unsigned long long  ulClockBefore, ulClockAfter;

		if (!isIntel() || !supportRDRAND())
		{
			// Generating error messages; we might want to fall back to LargerRand() or some
			//     other pseudo random number generator.
			fprintf(stderr, "The hardware does not make use of RDRAND instructions.\n");

			return 1;
		}

		uNumGoodSamples = 0;
		for (i = 0; i < NUM_RDRAND_ITERATIONS; i++) {
			uRetryCount = 0;

			ulClockBefore = __rdtsc();
			while (uRetryCount < RDRAND_RETRY)
			{
				if (_rdrand32_step(&uRdSample))
					break;
				uRetryCount++;
			}
			ulClockAfter = __rdtsc();

			// If we retried RDRAND_RETRY times without generating a good uRdSample, we
			//     skip over the statistics for this iteration.
			if (uRetryCount == RDRAND_RETRY)
				continue;

			// fRdSample is the uniformly distributed random number in [0, 1).
			fRdSample = (double)uRdSample/UINT_MAX;

			// Updating all the related statistics.
			fSumSamples += fRdSample;
			fSumSampleSquares += fRdSample*fRdSample;

			unsigned long   uElapsedCycles = (unsigned long) (ulClockAfter - ulClockBefore);
			fSumCycles += (double)uElapsedCycles;
			if (uElapsedCycles == uMinRDCycles)
			{
				uMinRDCount++;
			}
			else if (uElapsedCycles == uMaxRDCycles)
			{
				uMaxRDCount++;
			}
			else if (uElapsedCycles < uMinRDCycles)
			{
				// We have a new minimum cycle count from RDRAND.
				uMinRDCycles = uElapsedCycles;
				uMinRDCount = 1;
			}
			else if (uElapsedCycles > uMaxRDCycles)
			{
				// We have a new maximum cycle count from RDRAND.
				uMaxRDCycles = uElapsedCycles;
				uMaxRDCount = 1;
			}

			fNumRetries += (double)uRetryCount;
			uMaxRetries = (uRetryCount > uMaxRetries) ? uRetryCount : uMaxRetries;
			uNumGoodSamples++;
		}

		// Showing all the statistics.
		if (uNumGoodSamples)
		{
			// E(X) = Sum(X)/N
			fSumSamples /= uNumGoodSamples;
			// Var(X) = Sum(X^2)/N - E(X)*E(X)
			fSumSampleSquares = fSumSampleSquares / uNumGoodSamples - fSumSamples * fSumSamples;

			fSumCycles /= uNumGoodSamples;
			fNumRetries /= uNumGoodSamples;
			fprintf(stdout, "Mean of all %u samples: %lf\n", uNumGoodSamples, fSumSamples);
			fprintf(stdout, "Variance of all %u samples: %lf\n", uNumGoodSamples, fSumSampleSquares);
			fprintf(stdout, "On average each 32-bit RDRAND call costs %lf cycles.\n", fSumCycles);
			fprintf(stdout, "Minimum cycles RDRAND took is %lu cycles, occurring %lu times.\n", uMinRDCycles, uMinRDCount);
			fprintf(stdout, "Maximum cycles RDRAND took is %lu cycles, occurring %lu times.\n", uMaxRDCycles, uMaxRDCount);
			fprintf(stdout, "On average we had %lf%% retries for RDRAND; at most we invoked %u retries for one RDRAND.\n", fNumRetries*100, uMaxRetries);
		}
		else
		{
			fprintf(stderr, "We did not generate any good random number samples in %u tries.\n", NUM_RDRAND_ITERATIONS);
		}
	}
#endif
	// In the sixth part, we generate many normally distributed random variables with SIMD and showing the performance statistics.
	{
		double            fRdSample;
		double            fSumSamples = 0.0, fSumSampleSquares = 0.0;
		double            fSumCycles = 0.0;
		unsigned int      uNumGoodSamples;
		unsigned int      uMinRDCycles = UINT_MAX, uMinRDCount = 0;
		unsigned int      uMaxRDCycles = 0, uMaxRDCount = 0;
		unsigned long long  ulClockBefore, ulClockAfter;
		double (*GaussianRand_ptr)() = NULL;
		unsigned int	  uGaussianBinCounts[NUM_GAUSSIAN_BINS];
		double			  dfGaussianBinGap = 2.0 * GAUSSIAN_BIN_RANGE / NUM_GAUSSIAN_BINS;

		if (isIntel())
		{
			fprintf(stdout, "The running platform is based on Intel processors.\n\n\n");
		}
		else if (isAMD())
		{
			fprintf(stdout, "The running platform is based on AMD processors.\n\n\n");
		}
		else
		{
			fprintf(stdout, "The running platform is based on processors of an unknown brand.\n");
			return 1;
		}

		uNumGoodSamples = 0;
		// Using GaussianRandVec() when both AVX and AVX2 are supported.
		GaussianRand_ptr = supportAVX512F() ? GaussianRandVec : GaussianRand;

		for (i = 0; i < NUM_GAUSSIAN_BINS; i++)
			uGaussianBinCounts[i] = 0;

		for (i = 0; i < NUM_RDRAND_ITERATIONS; i++) {

			ulClockBefore = __rdtsc();
			fRdSample = GaussianRand_ptr();
			ulClockAfter = __rdtsc();

			// Updating all the related statistics.
			fSumSamples += fRdSample;
			fSumSampleSquares += fRdSample * fRdSample;

			// Which bin will dfGaussianRandNum fall in?
			int		iBinIndex = (int)((fRdSample + GAUSSIAN_BIN_RANGE) / dfGaussianBinGap);
			if ((iBinIndex >= 0) && (iBinIndex < NUM_GAUSSIAN_BINS))
				uGaussianBinCounts[iBinIndex]++;

			unsigned long   uElapsedCycles = (unsigned long)(ulClockAfter - ulClockBefore);
			fSumCycles += (double)uElapsedCycles;
			if (uElapsedCycles == uMinRDCycles)
			{
				uMinRDCount++;
			}
			else if (uElapsedCycles == uMaxRDCycles)
			{
				uMaxRDCount++;
			}
			else if (uElapsedCycles < uMinRDCycles)
			{
				// We have a new minimum cycle count from RDRAND.
				uMinRDCycles = uElapsedCycles;
				uMinRDCount = 1;
			}
			else if (uElapsedCycles > uMaxRDCycles)
			{
				// We have a new maximum cycle count from RDRAND.
				uMaxRDCycles = uElapsedCycles;
				uMaxRDCount = 1;
			}

			uNumGoodSamples++;
		}

		// Showing all the statistics.
		if (uNumGoodSamples)
		{
			// E(X) = Sum(X)/N
			fSumSamples /= uNumGoodSamples;
			// Var(X) = Sum(X^2)/N - E(X)*E(X)
			fSumSampleSquares = fSumSampleSquares / uNumGoodSamples - fSumSamples * fSumSamples;

			fSumCycles /= uNumGoodSamples;
			fprintf(stdout, "Mean of all %u samples: %lf\n", uNumGoodSamples, fSumSamples);
			fprintf(stdout, "Variance of all %u samples: %lf\n", uNumGoodSamples, fSumSampleSquares);
			fprintf(stdout, "On average each GaussianRandVec() call costs %lf cycles.\n", fSumCycles);
			fprintf(stdout, "Minimum cycles GaussianRandVec() took is %lu cycles, occurring %lu times.\n", uMinRDCycles, uMinRDCount);
			fprintf(stdout, "Maximum cycles GaussianRandVec() took is %lu cycles, occurring %lu times.\n", uMaxRDCycles, uMaxRDCount);
		}
		else
		{
			fprintf(stderr, "We did not generate any good random number samples in %u tries.\n", NUM_RDRAND_ITERATIONS);
		}
		fprintf(stdout, "\n\n\n");

		// Printing out the Gaussian bin counts.
		for (i = 0; i < NUM_GAUSSIAN_BINS; i++)
			fprintf(stdout, "%2dth bin (from %+f to %+f) contains %5d samples.\n", i,
				(-GAUSSIAN_BIN_RANGE + i * dfGaussianBinGap), (-GAUSSIAN_BIN_RANGE + (i + 1) * dfGaussianBinGap),
				uGaussianBinCounts[i]);
	}

	return 0;
}



