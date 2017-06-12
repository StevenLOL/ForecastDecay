// FUN.cpp: implementation of the FUN class.
//
//////////////////////////////////////////////////////////////////////

#include "FUN.h"
#include <iostream>
#include <math.h>
using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FUN::FUN()
{

}

FUN::~FUN()
{

}

//==========================计算参数v_excute程序==================================//
double *FUN::vParameter(double *ph, double *phn, int N, double f, double lined, double ae)
{
    double *vo;
	vo = new double[N - 2];
    double dan, dnb, hi1, hi2, d1, d2;
    int i;

    for (i = 1; i < N - 1; i++)
    {
        dan = lined / (N - 1) * i; //%到发射点的水平距离
        dnb = lined - dan;     // %到接受点的水平距离
        hi1 = fabs(ph[i] - ph[0]);
        hi2 = fabs(ph[i] - ph[N - 1]);
        d1 = sqrt(hi1 * hi1 + dan * dan);//%障碍物顶点到发射点的距离
        d2 = sqrt(hi2 * hi2 + dnb * dnb);//%障碍物顶点到接受的距离
        vo[i - 1] = 0.0316 * (phn[i - 1]) * sqrt(2 * (d1 + d2) / (300.0 / f * d1 * d2) * 1000);
    }
    return vo;

}
//================================================================================//
//==================fei_nie_R函数作用求采样点处的菲涅尔半径=======================//

double FUN::fei_nie_R(double *ph, int N, double f, double lined, int n)
{
    double d, k, a, d1, d2;
    double R;
    d = sqrt((ph[N - 1] - ph[0]) * (ph[N - 1] - ph[0]) + lined * lined) / 1000;//发射点到接受点的直线距离,单位是km
    k = fabs(ph[N - 1] - ph[0]) / lined;

    a = atan(k);               //倾角a
    d1 = (n - 1) * lined / (N - 1) / 1000 / cos(a);//%单位是km.
    d2 = d - d1;                // 单位是km.
    R = 550 * sqrt(d1 * d2 / d / f);  //样本点n处的菲涅尔半径
    return R;

}
//===============================================================================//
double QMax(double x,double y)
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

//=======================障碍物参数的确定，位置，最高点的菲涅尔余隙，宽度=====================//
void FUN::ObstacleParameter(double phn[], int N)
{
    //phn[]数组为菲涅尔余隙值，N为采样数列的总的点数。
    int i;

    //===========获得菲涅尔余隙中最大值zmax==========//

    zmax = phn[0];
    for (i = 1; i < N - 2; i++)
    {
        zmax = QMax(zmax, phn[i]);
    }
    //=======获得菲涅尔余隙最大值所在的位置nmax=======//
    for (i = 0; i < N - 2; i++)
    {
        if (phn[i] == zmax)
        {
            nmax = i + 1;
            break;
        }

    }
    //================计算障碍物的宽度===============//
    int m1 = 0, m2 = 0;
	if (zmax > 0)// 菲涅尔余隙大于0
	{
		if (nmax > 1)
		{
			for (i = nmax - 1; i >= 1; i--)
			{
				if (phn[i - 1] <= 0.2 * zmax)
				{
					m1 = i;
					break;
				}
				else
					m1 = i;
			}

			m1 = nmax - m1;
		}
		else
			m1 = 1;
		if (nmax < N - 2)
		{
			for (i = nmax + 1; i <= N - 2; i++)
			{
				if (phn[i - 1] <= 0.2 * zmax)
				{
					m2 = i;
					break;
				}
				else
					m2 = i;
			}

			m2 = m2 - nmax;
		}
		else
			m2 = 1;
	}

	else//菲涅尔余隙小于0
	{
		if (nmax > 1)
		{
			for (i = nmax - 1; i >= 1; i--)
			{
				if (phn[i - 1] <= 1.2 * zmax)
				{
					m1 = i;
					break;
				}
				else
					m1 = i;
			}

			m1 = nmax - m1;
		}
		else
			m1 = 0;
		if (nmax < N - 2)
		{
			for (i = nmax + 1; i <= N - 2; i++)
			{
				if (phn[i - 1] <= 1.2 * zmax)
				{
					m2 = i;
					break;
				}
				else
					m2 = i;
			}

			m2 = m2 - nmax;
		}
		else
			m2 = 0;
	}

	if (m1 == 0 || m2 == 0)
	{
		m = 0;
	}
	else
		m = QMax(m1, m2);

}
//=================================================================================//