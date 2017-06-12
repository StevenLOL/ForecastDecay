// Obstacle.cpp: implementation of the Obstacle class.
//
//////////////////////////////////////////////////////////////////////

#include "Obstacle.h"
#include "FUN.h"
#include <iostream>
#include <math.h>
using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Obstacle::Obstacle()
{

}

Obstacle::~Obstacle()
{

}
//===============================================================================//
double Omax(double x,double y)
{
	double result;
	if(x > y)
	{
		result = x;
	}
	else
	{
		result = y;
	}
	return result;
}
//===============================================================================//
void Obstacle::ObstacleNum(double phn[], double ph[], double lined, int N, double f)
{

    
    double R1, R2, R3;
    int p, q;
    int i;
	y5 = 0;
    R1 = R2 = R3 = 0.0;
    FUN Obstacle2;
    Obstacle2.ObstacleParameter(phn, N);
	Nmax[0] = 0;
	M[0] = 0;
	Nmax[2] = 0;
	M[2] = 0;
    Nmax[1] = Obstacle2.nmax;
    M[1] = Obstacle2.m;
   
    FUN FeiNieR2;
    R2 = FeiNieR2.fei_nie_R(ph, N, f, lined, Nmax[1] + 1);

    if (Obstacle2.zmax > (-0.6 * R2)) //障碍物的判断，障碍物处菲涅而余隙满足的条件。
    {
        y5 = y5 + 1;
    }
    p = Obstacle2.nmax - Obstacle2.m;
    q = Obstacle2.nmax + Obstacle2.m;

    
    
    int n1;
    if (p > 3)//满足此条件才可能&&前一段&&再出现满足条件的障碍物
    {
        double *phn1;
		phn1 = new double[Obstacle2.nmax - Obstacle2.m - 1];
        for (i = 0; i < Obstacle2.nmax - Obstacle2.m - 1; i++)
        {
            phn1[i] = phn[i];
        }

        n1 = Obstacle2.nmax - Obstacle2.m - 1; //前半段的余隙序列的点数
        FUN Obstacle1;
        Obstacle1.ObstacleParameter(phn1, n1 + 2); //  n1+2是总的采样点数
        FUN FeiNieR1;
        R1 = FeiNieR1.fei_nie_R(ph, N, f, lined, Obstacle1.nmax + 1);
        if (Obstacle1.zmax > (-0.6 * R1)) //障碍物的判断，障碍物处菲涅而余隙满足的条件。
		{
			if (Obstacle1.m>0)
			{
				y5 = y5 + 1;
			}

		}

        if (y5 == 1)
        {
            Obstacle1.nmax = 0;
            Obstacle1.m = 0;
        }
        Nmax[0] = Obstacle1.nmax;
        M[0] = Obstacle1.m;
    }
   

    if (q < N - 5)//满足此条件才可能&&前一段&&再出现满足条件的障碍物
    {
        double *phn2;
		phn2 = new double[N - 2 - Obstacle2.nmax - Obstacle2.m];
        for (i = 0; i < N - 2 - Obstacle2.nmax - Obstacle2.m; i++)
        {
            phn2[i] = phn[i + Obstacle2.nmax + Obstacle2.m];
        }
        FUN Obstacle3;
        Obstacle3.ObstacleParameter(phn2, N - 2 - Nmax[1] - M[1] + 2); //  n1+2是总的采样点数
        FUN FeiNieR3 ;
        R3 = FeiNieR3.fei_nie_R(ph, N, f, lined, Obstacle3.nmax + 1);
        Obstacle3.nmax = Obstacle3.nmax + q;


		if (Obstacle3.zmax > (-0.6 * R3)) //障碍物的判断，障碍物处菲涅而余隙满足的条件。
		{
			if (Obstacle3.m>0)
			{
				y5 = y5 + 1;
			}
		}

        if (y5 == 2 && Nmax[0] != 0)
        {
            Obstacle3.nmax = 0;
            Obstacle3.m = 0;
        }
        Nmax[2] = Obstacle3.nmax;

        M[2] = Obstacle3.m;
    }
   
    
}

//======================刀刃和圆形障碍物的判断，u>=3为刀刃，否则为圆形====================//
void Obstacle::U_Distinguish(double py[], double lined, int pointnum, double f, double ht, double hr, double ae)
{
    double dan, dnb;
    double *yn;
	yn = new double[pointnum];
    int i;
    //==============获得菲涅尔余隙值================//
    
    for (i = 1; i < pointnum + 1; i++)
    {
        dan = lined / (pointnum + 1) * i; //到发射点的水平距离
        dnb = lined - dan;      //到接受点的水平距离
        
        yn[i - 1] = py[i] + (dan * dnb / (2 * ae)) - (((py[0] + ht) * dnb + (py[pointnum + 1] + hr) * dan) / lined);
    }

    double H0, hh, hM, hM1, hN, hN1;
    int ymax, M, N;
    ymax = 1;
    M = 0; N = 0;
    //===========获得菲涅尔余隙中最大值ynmax==========//
    ynmax = yn[0];
    for (i = 1; i < pointnum; i++)
    {
        ynmax = Omax(ynmax, yn[i]);
    }
    //=======获得菲涅尔余隙最大值所在的位置nmax=======//
    for (i = 0; i < pointnum; i++)
    {
        if (yn[i] == ynmax)
        {
            ymax = i + 1;
            break;
        }

    }
    //==========获得位置nmax的菲涅尔半径===========//
    FUN FeiNie;
    Rm = FeiNie.fei_nie_R(py, pointnum + 2, f, lined, ymax + 1);

    H0 = Rm / (sqrt(3.0));

    hh = py[ymax] - H0;

    for (i = ymax; i > 0; i--)
    {
        if (py[i - 1] <= hh)
        {
            M = i;
           
            hM = py[M - 1] - hh;
            hM1 = py[M] - hh;
            if (fabs(hM) > fabs(hM1))
            {
                M = M + 1;
            }

            break;
        }
        else
            M = 1;

    }

    for (i = ymax + 2; i <= pointnum + 2; i++)
    {
        if (py[i - 1] <= hh)
        {
            N = i;
          
            hN1 = py[N - 2] - hh;
            hN = py[N - 1] - hh;
            if (fabs(hN) > fabs(hN1))
            {
                N = N - 1;
            }

            break;
        }
        else
            N = pointnum + 2;
    }
    //===================计算u值======================//
    double d, k, a, ro, X;
    
    d = sqrt((py[pointnum + 1] - py[0]) * (py[pointnum + 1] - py[0]) + lined * lined) / 1000;//发射点到接受点的直线距离,单位是km
    k = fabs(py[pointnum + 1] - py[0]) / lined;
    a = atan(k);

    ro = lined / (pointnum + 1) * (N - M) / cos(a);
    X = ymax * (lined / (pointnum + 1)) / lined;      //水平距离相除

    u = 2.02 * pow(X * (1 - X) * lined / ro, 2.0 / 3.0);
 
}