// attenuation_ra.h: interface for the attenuation_ra class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ATTENUATION_RA_H__1458CFF3_C531_4883_978C_199E9254EF9A__INCLUDED_)
#define AFX_ATTENUATION_RA_H__1458CFF3_C531_4883_978C_199E9254EF9A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class attenuation_ra  
{
public:
	attenuation_ra();
	virtual ~attenuation_ra();
	void shiju(double,double,double, double, double, double, double);
	void raosheParam(int, double [], double [], double, double, double);
	void raoshe(double, double, double, double, double, int, double[], double[], double, double, double);
	double atmoduct(double, double, double, double, double, double, double, double, double, double, double, double);
	void paramcalcu(double, double, double, double [], double [], double, int, double, double, double, double);
	double Radio_C(double**, double, double, double, double, double, int, int,double, double, double, double,double,double,double,double,double);

};

#endif // !defined(AFX_ATTENUATION_RA_H__1458CFF3_C531_4883_978C_199E9254EF9A__INCLUDED_)
