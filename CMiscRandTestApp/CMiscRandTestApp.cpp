// CMiscRandTestApp.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <stdlib.h>
#include "MiscRand.h"

#define	RAND_SEED			0
#define RAND_SEQUENCE_LEN	8
#define NUM_GAUSSIAN_RANDOM		524288
#define GAUSSIAN_BIN_RANGE		4			// We try to put the generated Gaussian random numbers into bins in
											// [-GAUSSIAN_BIN_RANGE, GAUSSIAN_BIN_RANGE]
#define NUM_GAUSSIAN_BINS		32

int main(int argc, char* argv[])
{
	unsigned long	i;
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
#endif
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
					(-GAUSSIAN_BIN_RANGE + i*dfGaussianBinGap), (-GAUSSIAN_BIN_RANGE + (i+1)*dfGaussianBinGap),
					uGaussianBinCounts[i]);
	}

	return 0;
}



