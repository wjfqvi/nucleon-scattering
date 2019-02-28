#include<complex>
#include<cmath>
#include "phase.h"
#include "functions.h"

using namespace std;


complex<double> Phase::mMatrix(int s,int m1,int m2,double acm)
{
	complex<double> mMatrix=complex<double>(0,0);
	complex<double> s11,s12,s22;
	complex<double> left1,left2,right1,right2;
	complex<double> mun,mcp;
	if(s==0){
		for(int j=0;j<=jMax;j++){
			mMatrix += (exp(2.0*iUnit*DToR*phaseun0[j])-1.0)*spherHms(j,0,cos(DToR*acm))*sqrt(2.0*j+1.0);
		}
		mMatrix = mMatrix*PiSqr/(iUnit*elabToKcm(elab));
		return mMatrix;
	}
	else if(s==1){

		mun = (exp(2.0*iUnit*DToR*phaseun1[0])-1.0)*spherHms(1,m2-m1,cos(DToR*acm))
					*cg_cff(1,1,0,m2-m1,m1)*cg_cff(1,1,0,0,m2)*sqrt(3.0);
		mcp = 0;

		mMatrix = mun;

		for(int j=1;j<=jMax;j++){

			left1=spherHms(j-1,m2-m1,cos(DToR*acm))*cg_cff(j-1,1,j,m2-m1,m1);
			left2=spherHms(j+1,m2-m1,cos(DToR*acm))*cg_cff(j+1,1,j,m2-m1,m1);

			s11=cos(DToR*2.0*e[j])*exp(DToR*2.0*iUnit*phasecpL[j])-1.0;
			s12=-iUnit*sin(DToR*2.0*e[j])*exp(DToR*iUnit*(phasecpL[j]+phasecpR[j]));
			s22=cos(DToR*2.0*e[j])*exp(DToR*2.0*iUnit*phasecpR[j])-1.0;

			right1=sqrt(2.0*(j-1)+1)*cg_cff(j-1,1,j,0,m2);
			right2=sqrt(2.0*(j+1)+1)*cg_cff(j+1,1,j,0,m2);

			mcp  = left1*(s11*right1+s12*right2);
			mcp += left2*(s12*right1+s22*right2);

			mun = (exp(DToR*2.0*iUnit*phaseun1[j])-1.0)*spherHms(j,m2-m1,cos(DToR*acm))*cg_cff(j,1,j,m2-m1,m1)
					*sqrt(2.0*j+1.0)*cg_cff(j,1,j,0,m2);

			mMatrix += mcp + mun;

		}

		mMatrix *= PiSqr/(iUnit*elabToKcm(elab));

		return mMatrix;
		
	}
}

