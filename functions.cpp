#include<cmath>


const double Pi=3.14;

double plgndr_s(int l,int m,double x)
{
	double plgndr_s;
	double pll,pmm,pmmp1,somx2;
	int ll,temp;
	pmm=1.0;
	if(m>=0 && l>=m && abs(x)<=1.0){
		if (m>0){
			somx2=sqrt((1.0-x)*(1.0+x));
			temp=1;
			for(int i=1;i<=m;i++){
				temp *= i*2-1;
			}
			pmm=temp*pow(somx2,m);
			if (m%2==1) pmm = -pmm;
		}
		if(l==m){plgndr_s=pmm;}
		else {
			pmmp1=x*(2*m+1)*pmm;
			if(l==m+1){
				plgndr_s=pmmp1;
			}
			else{
				for(int i=m+2;i<=l;i++){
					pll=(x*(2*i-1)*pmmp1-(i+m-1)*pmm)/(i-m);
						pmm=pmmp1;
						pmmp1=pll;
				}
				plgndr_s=pll;				
			} 
		}
	}
	else if(l<=abs(m) && abs(x)>1.0){plgndr_s=0;}
	return plgndr_s;
}


double spherHms(int l,int m,double x)
{

	double fact[30]={1,1,2,6,24,120,720,5040,40320,362880,3.6288e+06,3.99168e+07,4.79002e+08,
					6.22702e+09,8.71783e+10,1.30767e+12,2.09228e+13,3.55687e+14,6.40237e+15,
					1.21645e+17,2.4329e+18,5.10909e+19,1.124e+21,2.5852e+22,6.20448e+23,
					1.55112e+25,4.03291e+26,1.08889e+28,3.04888e+29,8.84176e+30};
	if(m>=0 && m<=l) {
		return sqrt((2*l+1.0)*fact[l-m]/(4*Pi*fact[l+m]));//*plgndr_s(l,m,x);
	}
	else if (m<0 && m>=-l){
		return pow((-1),abs(m))*sqrt((2*l+1)*fact[l-abs(m)]/((4*Pi)*fact[l+abs(m)]));//*plgndr_s(l,abs(m),x);
	}
}

int diracFunc(int l, int ls)
{
	if(l==ls) return 1;
	else return 0;
}


double cg_cff(int j1, int j2, int j3, int m1, int m2)
{
	double cg_cff,sumk,term;
	double fact[30]={1,1,2,6,24,120,720,5040,40320,362880,3.6288e+06,3.99168e+07,4.79002e+08,
					6.22702e+09,8.71783e+10,1.30767e+12,2.09228e+13,3.55687e+14,6.40237e+15,
					1.21645e+17,2.4329e+18,5.10909e+19,1.124e+21,2.5852e+22,6.20448e+23,
					1.55112e+25,4.03291e+26,1.08889e+28,3.04888e+29,8.84176e+30};
	int k,m3;
	m3=m1+m2;
	if(j3<abs(j1-j2) || j3>(j1+j2) || abs(m1)>j1 || abs(m2)>j2 || m3>j3 ){
		cg_cff = 0.0;
	}
	else{
		cg_cff=sqrt((j3*2+1)/fact[j1+j2+j3+1]);
		cg_cff=cg_cff*sqrt(fact[j1+j2-j3]*fact[j2+j3-j1]*fact[j3+j1-j2]);
		cg_cff=cg_cff*sqrt(fact[j1+m1]*fact[j1-m1]*fact[j2+m2]*fact[j2-m2]*fact[j3+m3]*fact[j3-m3]);
		sumk = 0.0;
		for(int i=0;i<26;i++){
			if((j1+j2-j3-i)>=0 && (j3-j1-m2+i)>=0 && (j3-j2+m1+i >=0) && (j1-m1-i)>=0 && (j2+m2-k)>=0 ){
				term=fact[j1+j2-j3-i]*fact[j3-j1-m2+i]*fact[j3-j2+m1+i]*fact[j1-m1-i]*fact[j2+m2-i]*fact[i];
				if(i%2==1) term=-term;
				sumk = sumk+1.0/term;
			}

		}
		cg_cff=cg_cff*sumk;
	}
}

