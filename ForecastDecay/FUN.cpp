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

//==========================�������v_excute����==================================//
double *FUN::vParameter(double *ph, double *phn, int N, double f, double lined, double ae)
{
    double *vo;
	vo = new double[N - 2];
    double dan, dnb, hi1, hi2, d1, d2;
    int i;

    for (i = 1; i < N - 1; i++)
    {
        dan = lined / (N - 1) * i; //%��������ˮƽ����
        dnb = lined - dan;     // %�����ܵ��ˮƽ����
        hi1 = fabs(ph[i] - ph[0]);
        hi2 = fabs(ph[i] - ph[N - 1]);
        d1 = sqrt(hi1 * hi1 + dan * dan);//%�ϰ��ﶥ�㵽�����ľ���
        d2 = sqrt(hi2 * hi2 + dnb * dnb);//%�ϰ��ﶥ�㵽���ܵľ���
        vo[i - 1] = 0.0316 * (phn[i - 1]) * sqrt(2 * (d1 + d2) / (300.0 / f * d1 * d2) * 1000);
    }
    return vo;

}
//================================================================================//
//==================fei_nie_R��������������㴦�ķ������뾶=======================//

double FUN::fei_nie_R(double *ph, int N, double f, double lined, int n)
{
    double d, k, a, d1, d2;
    double R;
    d = sqrt((ph[N - 1] - ph[0]) * (ph[N - 1] - ph[0]) + lined * lined) / 1000;//����㵽���ܵ��ֱ�߾���,��λ��km
    k = fabs(ph[N - 1] - ph[0]) / lined;

    a = atan(k);               //���a
    d1 = (n - 1) * lined / (N - 1) / 1000 / cos(a);//%��λ��km.
    d2 = d - d1;                // ��λ��km.
    R = 550 * sqrt(d1 * d2 / d / f);  //������n���ķ������뾶
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

//=======================�ϰ��������ȷ����λ�ã���ߵ�ķ�������϶�����=====================//
void FUN::ObstacleParameter(double phn[], int N)
{
    //phn[]����Ϊ��������϶ֵ��NΪ�������е��ܵĵ�����
    int i;

    //===========��÷�������϶�����ֵzmax==========//

    zmax = phn[0];
    for (i = 1; i < N - 2; i++)
    {
        zmax = QMax(zmax, phn[i]);
    }
    //=======��÷�������϶���ֵ���ڵ�λ��nmax=======//
    for (i = 0; i < N - 2; i++)
    {
        if (phn[i] == zmax)
        {
            nmax = i + 1;
            break;
        }

    }
    //================�����ϰ���Ŀ��===============//
    int m1 = 0, m2 = 0;
	if (zmax > 0)// ��������϶����0
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

	else//��������϶С��0
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