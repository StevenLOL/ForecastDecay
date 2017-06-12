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
	double u;    //障碍物类型参数u>=3为刀刃型，u<3为圆形
    double Rm;   //路径中采样点最大菲涅尔半径
    double ynmax;//精度更高的路径采样点中最大的菲涅尔余隙
	//==================障碍物个数和参数===========================//
    int y5;    //障碍物数量参数
    int Nmax[3]; //障碍物菲涅尔余隙最大位置参数数组
    int M[3];    //障碍物宽度参数数组
	void ObstacleNum(double [], double [], double, int, double);
	void U_Distinguish(double [], double, int, double, double, double, double);

};

#endif // !defined(AFX_OBSTACLE_H__03A4EE03_B83D_4374_8B10_C83E84EF95DF__INCLUDED_)
