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
	double P();
	double D();
	double AYY();
	double DT();
	double AP();
	double APT();
	double AT();
	double AXX();
	double AZX();
	double AZZ();
	double D0SK();
	double NNKK();
	double NNSK();
	double NSKN();
	double NSSN();
	double R();
	double RP();
	double RPT();
	double RT();
	double SGT();
};