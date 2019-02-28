/*
phaseun0-phasecpL-phaseun1-phasecpR-e
1S0     3P0
1P1 3S1 3P1 3D1 E1
1D2 3P2 3D2 3F2 E2
1F3 3D3 3F3 3G3 E3
1G4 3F4 3G4 3H4 E4
1H5 3G5 3H5 3I5 E5
...
*/

#ifndef PHASE_H_
#define PHASE_H_

#include<complex>
#include<cmath>
#include "constant.h"

using namespace std;

class Phase{
// private:
public:
	int jMax;
	double phaseun0[10];
	double phaseun1[10];
	double phasecpL[10];
	double phasecpR[10];
	double e[10];
	double acm;
	double elab;
// public:
	complex<double> mMatrix(int s,int m,int mp,double acm);

};

#endif
