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
		//===================�����ϰ������=====================================//
    int nmax;     //�ϰ����������϶����λ��
    double zmax;  //�ϰ������ķ�������϶ֵ
    int m;        //�ϰ���Ŀ��
	double *vParameter(double[], double[], int, double, double,double);
	double fei_nie_R(double[], int, double, double, int);
	void ObstacleParameter(double[], int);

};

#endif // !defined(AFX_FUN_H__88B961F1_D739_43C0_A99F_0762057B9AC4__INCLUDED_)
