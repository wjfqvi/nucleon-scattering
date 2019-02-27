#ifndef SACLAY_H_
#define SACLAY_H_

#include<complex>

using namespace std;

class Saclay{
// private:

public:
	complex<double> a,b,c,d,e;
	double acm,alpha,beta;
	const double DToR=3.14159/180.0;

	double I0000();
	double DSG();
	
	double AYY();
	double D();
	double DT();
	double P();

	double R();
	double RP();	
	double AP();

	double RT();
	double AT();
	double RPT();
	double APT();

	double D0SK();

	double AXX();
	double AZX();
	double AZZ();

	double NSKN();
	double NSSN();
	double NNKK();
	double NNSK();

	double SGT();
};

#endif
