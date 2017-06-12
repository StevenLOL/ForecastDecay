// FUN.h: interface for the FUN class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FUN_H__88B961F1_D739_43C0_A99F_0762057B9AC4__INCLUDED_)
#define AFX_FUN_H__88B961F1_D739_43C0_A99F_0762057B9AC4__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class FUN  
{
public:
	FUN();
	virtual ~FUN();
		//===================单个障碍物参数=====================================//
    int nmax;     //障碍物菲涅尔余隙最大的位置
    double zmax;  //障碍物最大的菲涅尔余隙值
    int m;        //障碍物的宽度
	double *vParameter(double[], double[], int, double, double,double);
	double fei_nie_R(double[], int, double, double, int);
	void ObstacleParameter(double[], int);

};

#endif // !defined(AFX_FUN_H__88B961F1_D739_43C0_A99F_0762057B9AC4__INCLUDED_)
