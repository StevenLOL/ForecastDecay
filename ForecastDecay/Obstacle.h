// Obstacle.h: interface for the Obstacle class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_OBSTACLE_H__03A4EE03_B83D_4374_8B10_C83E84EF95DF__INCLUDED_)
#define AFX_OBSTACLE_H__03A4EE03_B83D_4374_8B10_C83E84EF95DF__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class Obstacle  
{
public:
	Obstacle();
	virtual ~Obstacle();
	double u;    //�ϰ������Ͳ���u>=3Ϊ�����ͣ�u<3ΪԲ��
    double Rm;   //·���в��������������뾶
    double ynmax;//���ȸ��ߵ�·�������������ķ�������϶
	//==================�ϰ�������Ͳ���===========================//
    int y5;    //�ϰ�����������
    int Nmax[3]; //�ϰ����������϶���λ�ò�������
    int M[3];    //�ϰ����Ȳ�������
	void ObstacleNum(double [], double [], double, int, double);
	void U_Distinguish(double [], double, int, double, double, double, double);

};

#endif // !defined(AFX_OBSTACLE_H__03A4EE03_B83D_4374_8B10_C83E84EF95DF__INCLUDED_)
