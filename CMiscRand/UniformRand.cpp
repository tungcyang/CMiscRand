#include "MiscRand.h"

static unsigned long		uLargerRandSeed = 0;

void	__cdecl	sLargerRand(unsigned long _Seed)
{
	uLargerRandSeed = _Seed;
}

long		__cdecl	LargerRand()
{
	// We no longer compute ((uLargerRandSeed = uLargerRandSeed * 214013L + 2531011L) >> 16) & 0x7fff;
	// we only evalute the linear congruential generator and the users need to apply the
	// right shift and the bitmask.
	return (long)(uLargerRandSeed = uLargerRandSeed * 214013L + 2531011L);
}