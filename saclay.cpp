#include "saclay.h"


double Saclay::I0000()
{
	return (norm(a)+norm(b)+norm(c)+norm(d)+norm(e))/2;
}

double Saclay::DSG()
{
	return I0000()*389387.0;
}

double Saclay::AYY()
{
	return (norm(a)-norm(b)-norm(c)+norm(d)+norm(e))/2/I0000();
}

double Saclay::D()
{
	return (norm(a)+norm(b)-norm(c)-norm(d)+norm(e))/2/I0000();
}

double Saclay::DT()
{
	return (norm(a)-norm(b)+norm(c)-norm(d)+norm(e))/2/I0000();
}

double Saclay::P()
{
	return real(conj(a)*e)/I0000();
}

double Saclay::R()
{
	return  (real(conj(a)*b)*cos(DToR*(alpha+acm/2))
			+real(conj(c)*d)*sin(DToR*(alpha-acm/2))
			-imag(conj(b)*e)*sin(DToR*(alpha-acm/2)))/I0000();
}

double Saclay::RP()
{
	return  (real(conj(a)*b)*sin(DToR*(alpha+acm/2))
			+real(conj(c)*d)*sin(DToR*(alpha-acm/2))
			+imag(conj(b)*e)*cos(DToR*(alpha+acm/2)))/I0000();
}

double Saclay::AP()
{
	return  (real(conj(a)*b)*cos(DToR*(alpha+acm/2))
			-real(conj(c)*d)*cos(DToR*(alpha-acm/2))
			-imag(conj(b)*e)*sin(DToR*(alpha+acm/2)))/I0000();
}

double Saclay::RT()
{
	return (-real(conj(a)*c)*cos(DToR*(beta+acm/2))
			-real(conj(b)*d)*cos(DToR*(beta-acm/2))
			-imag(conj(c)*e)*cos(DToR*(beta+acm/2)))/I0000();
}

double Saclay::AT()
{
	return (-real(conj(a)*c)*sin(DToR*(beta+acm/2))
			+real(conj(b)*d)*sin(DToR*(beta-acm/2))
			+imag(conj(c)*e)*sin(DToR*(beta+acm/2)))/I0000();
}

double Saclay::RPT()
{
	return (-real(conj(a)*c)*sin(DToR*(beta+acm/2))
			-real(conj(b)*d)*sin(DToR*(beta-acm/2))
			-imag(conj(c)*e)*cos(DToR*(beta+acm/2)))/I0000();
}

double Saclay::APT()
{
	return  (real(conj(a)*c)*cos(DToR*(beta+acm/2))
			-real(conj(b)*d)*cos(DToR*(beta-acm/2))
			+imag(conj(b)*e)*sin(DToR*(beta+acm/2)))/I0000();
}

double Saclay::D0SK()
{
	return ( real(conj(a)*b)*sin(DToR*(beta+acm/2))
			-real(conj(c)*d)*sin(DToR*(beta-acm/2))
			+imag(conj(b)*e)*cos(DToR*(beta+acm/2)))/I0000();
}

double Saclay::AXX()
{
	return real(conj(a)*d)*cos(DToR*acm)+real(conj(b)*c)-imag(conj(d)*e)*sin(DToR*acm);
}

double Saclay::AZX()
{
	return real(conj(a)*d)*sin(DToR*acm)+imag(conj(d)*e)*cos(DToR*acm);
}

double Saclay::AZZ()
{
	return -real(conj(a)*d)*cos(DToR*acm)+real(conj(b)*c)+imag(conj(d)*e)*sin(DToR*acm);
}

double Saclay::NSKN()
{
	return ( real(conj(c)*e)*sin(DToR*(beta+acm/2))
			-imag(conj(a)*e)*cos(DToR*(beta+acm/2))
			+imag(conj(b)*d)*cos(DToR*(beta-acm/2)))/I0000();
}

double Saclay::NSSN()
{
	return (-real(conj(c)*e)*cos(DToR*(beta+acm/2))
			-imag(conj(a)*c)*sin(DToR*(beta+acm/2))
			-imag(conj(b)*d)*sin(DToR*(beta-acm/2)))/I0000();
}

double Saclay::NNKK()
{
	return -real(conj(d)*e)*cos(DToR*acm)-imag(conj(a)*d)*sin(DToR*acm);
}

double Saclay::NNSK()
{
	return -real(conj(d)*e)*sin(DToR*acm)+imag(conj(a)*d)*cos(DToR*acm)+imag(conj(b)*c);
}







