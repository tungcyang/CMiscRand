#pragma once
// isIntel() returns 1 if the code is executed on an Intel CPU, 0 otherwise.
bool	isIntel();

// isAMD() returns 1 if the code is executed on an AMD CPU, 0 otherwise.
bool	isAMD();

// supportRDRAND() returns TRUE if the CPU where the code is executed supports
// RDRAND feature.
bool	supportRDRAND();

// supportRDSEED() returns TRUE if the CPU where the code is executed supports
// RDSEED feature.
bool	supportRDSEED();

// supportAVX() returns TRUE if the CPU where the code is executed supports
// AVX feature.
bool	supportAVX();

// supportAVX2() returns TRUE if the CPU where the code is executed supports
// AVX2 feature.
bool	supportAVX2();