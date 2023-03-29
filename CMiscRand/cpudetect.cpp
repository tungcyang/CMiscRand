#include <stdio.h>
#include "intrin.h"
#include "cpudetect.h"


// isIntel() returns 1 if the code is executed on an Intel CPU, 0 otherwise.  Please refer to
//     https://en.wikipedia.org/wiki/CPUID
// for more information.
bool	isIntel()
{
	const int	Genu = 0x756E6547;  // "uneG", little-endian for [Genu]ineIntel.
	const int   Inei = 0x49656E69;  // "Ieni", little-endian for Genu[ineI]ntel.
	const int	ntel = 0x6C65746E;  // "letn", little-endian for GenuineI[ntel].

	int		cpuInfo[4];

	__cpuid(cpuInfo, 0);
	// EBX, EDX and ECX are returned in cpuInfo[1-3] when input EAX is 0.
	if ((cpuInfo[1] == Genu) && (cpuInfo[2] == ntel) && (cpuInfo[3] == Inei))
		return 1;
	else
		return 0;
}

// isAMD() returns 1 if the code is executed on an AMD CPU, 0 otherwise.  Please refer to
//     https://en.wikipedia.org/wiki/CPUID
// for more information.
bool	isAMD()
{
	const int	Auth = 0x68747541;  // "htuA", little-endian for [Auth]enticAMD.
	const int   enti = 0x69746E65;  // "itne", little-endian for Auth[enti]cAMD.
	const int	cAMD = 0x444D4163;  // "DMAc", little-endian for Authenti[cAMD].

	int		cpuInfo[4];

	__cpuid(cpuInfo, 0);
	// EBX, EDX and ECX are returned in cpuInfo[1-3] when input EAX is 0.
	if ((cpuInfo[1] == Auth) && (cpuInfo[2] == cAMD) && (cpuInfo[3] == enti))
		return 1;
	else
		return 0;
}

// supportRDRAND() returns TRUE if the CPU where the code is executed supports
// RDRAND feature.  Please refer to
//     https://en.wikipedia.org/wiki/CPUID
// for more information.
bool	supportRDRAND()
{
	int		cpuInfo[4];

	__cpuid(cpuInfo, 1);
	// EBX, ECX and EDX are returned in cpuInfo[1-3].
	return !!(cpuInfo[2] & 0x40000000);
}

// supportRDSEED() returns TRUE if the CPU where the code is executed supports
// RDSEED feature.  Please refer to
//     https://en.wikipedia.org/wiki/CPUID
// for more information.
bool	supportRDSEED()
{
	int		cpuInfo[4];

	__cpuid(cpuInfo, 7);
	// EBX, ECX and EDX are returned in cpuInfo[1-3].
	return !!(cpuInfo[1] & 0x00040000);
}

// supportAVX() returns TRUE if the CPU where the code is executed supports
// AVX feature.  Please refer to
//     https://en.wikipedia.org/wiki/CPUID
// for more information.
bool	supportAVX()
{
	int		cpuInfo[4];

	__cpuid(cpuInfo, 1);
	// EBX, ECX and EDX are returned in cpuInfo[1-3].
	return !!(cpuInfo[2] & 0x10000000);
}

// supportAVX2() returns TRUE if the CPU where the code is executed supports
// AVX2 feature.  Please refer to
//     https://en.wikipedia.org/wiki/CPUID
// for more information.
bool	supportAVX2()
{
	int		cpuInfo[4];

	__cpuid(cpuInfo, 7);
	// EBX, ECX and EDX are returned in cpuInfo[1-3].
	return !!(cpuInfo[1] & 0x00000020);
}