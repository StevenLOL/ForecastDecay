// ITUrp.cpp: implementation of the ITUrp class.
//
//////////////////////////////////////////////////////////////////////

#include "ITUrp.h"
#include <iostream>
#include <math.h>
#include "Obstacle.h"
#include "FUN.h"
using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ITUrp::ITUrp()
{

}

ITUrp::~ITUrp()
{

}

 double OutVisual(double ph[], int LinePoint, double Distance, double fre, double T_high, double R_high, int y3, int y4, double aa, int ee, double aee)
{
    //OutVisual函数的功能是计算平滑地球视距外的绕射衰减
    //ph数组是路径采样点的海拔高度值， LinePoint是路径总的采样点数， Distance是发射机到接收机的距离，fre是信号频率，T_high是发射机的天线高度，R_high是接收机的天线高度
    //y3参数是极化参数，y3=1代表水平极化，y3=0代表垂直极化，y4参数是地形参数，y4=0代表陆地，y4=1代表海洋，aa是有效电导率，ee是相对介电常数，aee是地球有效半径。
    double k = 4.0 / 3.0;//这个数据和天气情况有关。4/3是最大的数
    double ae1 = (double)aee / 1000;//%单位是km
    double Kh, K, B, d, ht1, hr1;

    double X, Y1, Y2, Fx, Gy1, Gy2;
    double Fr;
    Kh = 0.36 * pow((ae1 * fre), (float)(-1.0 / 3.0)) * pow((float)((ee - 1) * (ee - 1) + (18000 * aa / fre) * (18000 * aa / fre)), (float)(-1.0 / 4.0));
    if (y3 == 1)
        K = Kh;
    else
        K = Kh * sqrt(ee * ee + (18000 * aa / fre) * (18000 * aa / fre));

    if (K < 0.001)
        Fr = 0;
    else
    {
        if (y3 == 1 || (y3 == 0 && (y4 == 0 && fre >= 20) || y4 == 1 && fre >= 300))
            B = 1;
        else
        {
            K = sqrt(6.89 * aa / pow(k, (float)(2.0 / 3.0)) / pow(fre, (float)(5.0 / 3.0)));
            B = (1 + 1.6 * K * K + 0.75 * K * K * K * K) / (1 + 4.5 * K * K + 1.35 * K * K * K * K);
        }
        
        ht1 = ph[0] + T_high;
        hr1 = ph[LinePoint - 1] + R_high;
        d = sqrt(Distance * Distance + (ht1 - hr1) * (ht1 - hr1)) / 1000;

        X = 2.2 * B * pow(fre, (float)(1.0 / 3.0)) * pow(ae1, (float)(-2.0 / 3.0)) * d;//%归一化长度
        Y1 = 9.6 * 0.001 * B * pow(fre, (float)(2.0 / 3.0)) * pow(ae1, (float)(-1.0 / 3.0)) * ph[0];
        Y2 = 9.6 * 0.001 * B * pow(fre, (float)(2.0 / 3.0)) * pow(ae1, (float)(-1.0 / 3.0)) * (ph[LinePoint - 1]);
        Fx = 11 + 10 * log10(X) - 17.6 * X;

        if (Y1 > 2)
            Gy1 = 17.6 * sqrt(Y1 - 1.1) - 5 * log10(Y1 - 1.1) - 8;
        else if (10 * K < Y1 && Y1 < 2)
            Gy1 = 20 * log10(Y1 + 0.1 * Y1 * Y1 * Y1);
        else if (K / 10 < Y1 && Y1 < 10 * K)
            Gy1 = 2 + 20 * log10(K) + 9 * log10(Y1 / K) * (log10(Y1 / K) + 1);
        else if (Y1 < K / 10)
            Gy1 = 2 + 20 * log10(K);
        else
            Gy1 = 0;


        if (Y2 > 2)
            Gy2 = 17.6 * sqrt(Y2 - 1.1) - 5 * log10(Y2 - 1.1) - 8;
        else if (10 * K < Y2 && Y2 < 2)
            Gy2 = 20 * log10(Y2 + 0.1 * Y2 * Y2 * Y2);
        else if (K / 10 < Y2 && Y2 < 10 * K)
            Gy2 = 2 + 20 * log10(K) + 9 * log10(Y2 / K) * (log10(Y2 / K) + 1);
        else if (Y2 < K / 10)
            Gy2 = 2 + 20 * log10(K);
        else
            Gy2 = 0;

        Fr = Fx + Gy1 + Gy2;
    }
    return Fr;
}
//=======================================================================================================================//
//===========================光滑地球视距内衰减损耗计算==============================//
double InVisual(double ph[], int N, double lined, double f, double ht, double hr, int y3, int y4, double a, int e, double ae)
{
    double ae1 = ae / 1000;//单位是km
  
    double h1, h2, hm, m, c, b, d1, d2, h, hreq;
    double A, Ah, aem, aem1;
    h1 = ph[0] + ht;
    h2 = ph[N - 1] + hr;
    hm = sqrt(h1) + sqrt(h2);
    m = lined * lined / (4 * ae * (h1 + h2));
    c = fabs(h1 - h2) / (h1 + h2);
    b = 2 * sqrt((m + 1) / m / 3) * cos(3.1416 / 3 + 1.0 / 3.0 * acos(1.5 * c * sqrt(3 * m / ((m + 1) * (m + 1) * (m + 1)))));
    d1 = lined / 2 * (1 + b);
    d2 = lined - d1;
    h = ((h1 - d1 * d1 / (2 * ae)) * d2 + (h2 - d2 * d2 / (2 * ae)) * d1) / lined;
    hreq = 0.552 * sqrt(d1 * d2 * 300.0 / f / lined);
    if (h > hreq)
        A = 0;
    else
    {
        aem = 0.5 * (lined / hm) * (lined / hm);//单位是km
        aem1 = aem * 1000;
        Ah = OutVisual(ph, N, lined, f, ht, hr, y3, y4, a, e, aem1);
        if (Ah < 0)
            A = 0;
        else
            A = (1 - h / hreq) * Ah;

    }
    return A;
}
//=======================================================================================================================//
//===============================================================================//
double max(double x,double y)
{
	double result=(x>y)?x:y;

	return result;
}
double min(double x,double y)
{
	double result=(x<y)?x:y;

	return result;
}
 //========================单个刀刃障碍物衰减损耗计算===============================//
double L_danfeng(double ph[], double phn[], int N, double lined, double f, double ae)
{
    double *v;
	v = new double[N - 2];
    double vmax, Jvp, L_dB;
    int i;
    FUN vPara;
    v = vPara.vParameter(ph, phn, N, f, lined, ae);   //获取参数v的首地址。
    
    vmax = v[0];
    for (i = 1; i < N - 2; i++)
    {
        vmax = max(vmax, v[i]);          //获取v中的最大值vmax
    }

    if (vmax > -0.78)
        Jvp = 6.9 + 20 * log10(sqrt((vmax - 0.1) * (vmax - 0.1) + 1) + vmax - 0.1);
    else
        Jvp = 0;

    L_dB = Jvp;

	delete[]v;
    return L_dB;


}

//=======================================================================================================================//
//========================双重刀刃障碍物衰减损耗计算===============================//
double L_doublefeng(double ph[], double phn[], int N, double lined, double f, double ae, int nmax[]) 
{
    double dd, a, b, c;
    int n1, n2;
    n1 = 0; n2 = 0;
    dd = lined / (N - 1);
    if (nmax[0] == 0)
    {
        n1 = nmax[1] + 1;
        n2 = nmax[2] + 1;
    }
    else if (nmax[2] == 0)
    {
        n1 = nmax[0] + 1;
        n2 = nmax[1] + 1;

    }
    a = (n1 - 1) * dd; b = (n2 - n1) * dd; c = (N - n2) * dd;

    /////////////////////第一个障碍物的衰减//////////////////////
    double h11, h12, d11, d12, h1, vo1, Jvp1, L1;
    
    h11 = fabs(ph[n1 - 1] - ph[0]);
    h12 = fabs(ph[n1 - 1] - ph[n2 - 1]);
    d11 = sqrt(h11 * h11 + a * a);      //障碍物顶点到发射点的距离
    d12 = sqrt(h12 * h12 + b * b);      //障碍物顶点到接受的距离
   
    h1 = ph[n1 - 1] + (a * b / (2 * ae)) - (ph[0] * b + ph[n2 - 1] * a) / (a + b);
    vo1 = 0.0316 * h1 * sqrt(2 * (d11 + d12) / (300.0 / f * d11 * d12) * 1000);
    if (vo1 > -0.78)
        Jvp1 = 6.9 + 20 * log10(sqrt((vo1 - 0.1) * (vo1 - 0.1) + 1) + vo1 - 0.1);
    else
        Jvp1 = 0;

    L1 = Jvp1;
    /////////////////////第二个障碍物的衰减//////////////////////
    double h21, h22, d21, d22, h2, vo2, Jvp2, L2;
    
    h21 = fabs(ph[n2 - 1] - ph[n1 - 1]);
    h22 = fabs(ph[N - 1] - ph[n2 - 1]);
    d21 = sqrt(h21 * h21 + b * b);//障碍物顶点到发射点的距离
    d22 = sqrt(h22 * h22 + c * c);//障碍物顶点到接受的距离
    
    h2 = ph[n2 - 1] + (b * c / (2 * ae)) - ((ph[n1 - 1] * c + ph[N - 1] * b) / (b + c));
    vo2 = 0.0316 * h2 * sqrt(2 * (d21 + d22) / (300.0 / f * d21 * d22) * 1000);
    if (vo2 > -0.78)
        Jvp2 = 6.9 + 20 * log10(sqrt((vo2 - 0.1) * (vo2 - 0.1) + 1) + vo2 - 0.1);
    else
        Jvp2 = 0;
    L2 = Jvp2;
    //////////////////////////////////////////////////////////
    double hh, L_dB, Lc, p, q, aa, Tc;
    double r1, r2, kk, h13, h31;

    double lamda = 300.0 / f;
    if (ph[n1 - 1] > ph[n2 - 1])
        hh = (ph[n1 - 1]) / (ph[n2 - 1]);
    else

        hh = (ph[n2 - 1]) / (ph[n1 - 1]);

    if (hh < 1.22) //两个障碍物相当
    {
        Lc = 10 * log10((a + b) * (b + c) / b / (a + b + c));
        if (L1 > 15 && L2 > 15)//当满足这个条件时，Lc修正项才有效
            L_dB = L1 + L2 + Lc;
        else
            L_dB = L1 + L2;
    }
    else  //主副障碍物衰减计算
    {

       
        p = sqrt(2 / lamda * (a + b + c) / (a * (b + c))) * (phn[n1 - 2]);
        q = sqrt(2 / lamda * (a + b + c) / (c * (b + a))) * (phn[n2 - 2]);
        aa = atan(sqrt(b * (a + b + c) / c / a));
        if (q / p < 0)
        { Tc = -8.0; }
        else
        { Tc = (12 - 20 * log10(2 / (1 - aa / 3.1416))) * pow(q / p, 2 * p); }
        // ========菲涅尔半径的计算==========//
        FUN FeiNieR1;
        r1 = FeiNieR1.fei_nie_R(ph, N, f, lined, n1);
        FUN FeiNieR2;
        r2 = FeiNieR2.fei_nie_R(ph, N, f, lined, n2);
       
        kk = (phn[n1 - 2]) / r1 - (phn[n2 - 2]) / r2;
        if (kk > 0)//%第一个障碍物为主障碍物
        {
            
            h13 = fabs(ph[n1 - 1] - ph[N - 1]);
            d11 = sqrt(h11 * h11 + a * a);//%障碍物顶点到发射点的距离
            d12 = sqrt(h13 * h13 + (b + c) * (b + c));//%障碍物顶点到接受的距离
           
            h1 = ph[n1 - 1] + (a * (b + c) / (2 * ae)) - ((ph[0] * (b + c) + (ph[N - 1] * a)) / (a + b + c));
            vo1 = 0.0316 * h1 * sqrt(2 * (d11 + d12) / (lamda * d11 * d12) * 1000);
            if (vo1 > -0.78)
                Jvp1 = 6.9 + 20 * log10(sqrt((vo1 - 0.1) * (vo1 - 0.1) + 1) + vo1 - 0.1);

            else
                Jvp1 = 0;

            L1 = Jvp1;
            L_dB = L1 + L2 - Tc;
        }
        //第二个障碍物为主障碍物时衰减
        else
        {
            
            h31 = fabs(ph[n2 - 1] - ph[0]);
            d21 = sqrt(h31 * h31 + (a + b) * (a + b));
            d22 = sqrt(h22 * h22 + c * c);
            
            h2 = ph[n2 - 1] + ((a + b) * c / (2 * ae)) - ((ph[0] * c + (ph[N - 1]) * (a + b)) / (a + b + c));
            vo2 = 0.0316 * h2 * sqrt(2 * (d21 + d22) / (lamda * d21 * d22) * 1000);
            if (vo2 > -0.78)
                Jvp2 = 6.9 + 20 * log10(sqrt((vo2 - 0.1) * (vo2 - 0.1) + 1) + vo2 - 0.1);
            else
                Jvp2 = 0;

            L2 = Jvp2;
            L_dB = L1 + L2 - Tc;
        }
        
        L_dB = fabs(L_dB);   
    }
    return L_dB;
}
//===============================================================================//
//========================三个刀刃障碍物衰减损耗计算===============================//
double L_duofeng(double ph[], double phn[], int N, double lined, double f, double ht, double hr, double ae)
{
    double *v; 
	v = new double[N - 2];
    double vmax, Jvp, Jvt, Jvr, L_dB;
    int i;
    FUN vPara;
    v = vPara.vParameter(ph, phn, N, f, lined, ae);//获得参数v

    //===============主障碍物的衰减损耗计算=================//
    vmax = v[0];
    for (i = 1; i < N - 2; i++)
    {
        vmax = max(vmax, v[i]);          //获取v中的最大值vmax
    }

    int maxp;
    maxp = 0;
    for (i = 0; i < N - 2; i++)
    {
        if (v[i] == vmax)
        {
            maxp = i + 1;                  //v值最大的位置
            break;
        }
    }


    if (vmax > -0.78)//绕射系数最大点>-0.78,计算主峰绕射损耗Jvp,分别对最大点（主峰）两端路径进行同样的分析
    {
        Jvp = 6.9 + 20 * log10(sqrt((vmax - 0.1) * (vmax - 0.1) + 1) + vmax - 0.1);

        //===============第一个障碍物的衰减损耗计算=================//


        double lined1 = lined / (N - 1) * maxp;
        double *v1;
		v1 = new double[maxp - 1]; //前一段的v的点数为maxp-1,对应的h序列为maxp+1;
        
        double *phnt;
		phnt=new double[maxp - 1];//第一个障碍物的菲涅尔余隙
       
        double dant, dnbt;
        for (i = 1; i < maxp; i++)
        {
            dant = lined / (N - 1) * i;           //到发射点的水平距离
            dnbt = lined1 - dant;              //到接受点的水平距离
            
            phnt[i - 1] = ph[i] + (dant * dnbt / (2 * ae)) - (((ph[0] + ht) * dnbt + (ph[maxp] + hr) * dant) / lined1);  //菲涅尔余隙计算
           
        }
        v1 = vPara.vParameter(ph, phnt, maxp + 1, f, lined1, ae);//获得参数v
       
        ///////////////////////////////////////////////////

        double vtmax;
        if (maxp - 1 < 4)
            vtmax = -10;
        else
        {
            vtmax = v1[0];

            for (i = 1; i < maxp - 1; i++)
            {
                vtmax = max(vtmax, v1[i]);          //获取前半段序列v中的最大值vtmax
            }
            
        }

        if (vtmax > -0.78)

            Jvt = 6.9 + 20 * log10(sqrt((vtmax - 0.1) * (vtmax - 0.1) + 1) + vtmax - 0.1);
        else
            Jvt = 0;

        //===============第二个障碍物的衰减损耗计算=================//
       
        if (N - maxp > 2)
        {
            double lined2 = lined / (N - 1) * (N - maxp - 1);
            double *v2;
			v2 = new double[N - 2 - maxp]; //后一段的v的点数为N-2-maxp,对应的h序列为N-maxp;
            double *phr;
			phr = new double[N - maxp];
            double *phnr;
			phnr = new double[N - 2 - maxp];
            double danr, dnbr;
            for (i = 0; i < (N - maxp); i++)
            {
                phr[i] = ph[i + maxp];

            }
            for (i = 1; i < N - maxp - 1; i++)
            {
                danr = lined / (N - 1) * i;           //到发射点的水平距离
                dnbr = lined2 - danr;              //到接受点的水平距离
                phnr[i - 1] = phr[i] + (danr * dnbr / (2 * ae)) - (((ph[maxp] + ht) * dnbr + (ph[N - 1] + hr) * danr) / lined2);  //菲涅尔余隙计算
               
            } 
            v2 = vPara.vParameter(phr, phnr, N - maxp, f, lined2, ae);//获得参数v

            double vrmax;
            if (N - 2 - maxp < 4)
            {
                vrmax = -10;
                
            }
            else
            {

                vrmax = v2[0];

                for (i = 1; i < N - 2 - maxp; i++)
                {
                    vrmax = max(vrmax, v2[i]);          //获取后半段序列v中的最大值vrmax
                }
                
            }
            if (vrmax > -0.78)
                Jvr = 6.9 + 20 * log10(sqrt((vrmax - 0.1) * (vrmax - 0.1) + 1) + vrmax - 0.1);
            else
                Jvr = 0;

			delete[]v2;
			delete[]phr;
			delete[]phnr;

        }
        else
            Jvr = 0;
        //====根据两个辅助峰和主峰绕射损耗及其位置进行修正得到总的绕射损耗====//


        double T, C;

        T = 1 - exp(-1 * Jvp / 6);

        C = 10 + 0.04 * lined / 1000;
        L_dB = Jvp + T * (Jvt + Jvr + C);

		delete[]v1;
		delete[]phnt;

    }
    else
        L_dB = 0;

	delete[]v;
    return L_dB;
}

 //===========================单个圆形障碍物的衰减计算==============================//
double yuan_u_L(double u, double Rm, double ynmax)
{
    double Hc, Vo, Ar, Ho, L_dB;
    Ar = 0.0;
    if (u <= 0.79)
        Ar = 5.5;
    else if (u > 0.79 && u <= 1.09)
        Ar = 3.3;
    else if (u > 1.09 && u <= 1.9)
        Ar = 2.0;
    else if (u > 1.9 && u <= 2.2)
        Ar = 1.8;
    else if (u > 2.2 && u < 3)
        Ar = 1.6;

    Hc = ynmax;
    Ho = -Rm / sqrt(3.0);//自由空间传播路径余隙  
    Vo = 14.42 + Ar * pow(fabs(u - 1.4), 1.5) - 20 * log10(u);
    L_dB = Vo * (1 - Hc / Ho);
    if (L_dB < 0)
        L_dB = 0.0;

    return L_dB;

}
//===============================================================================//
//==========================多圆障碍物子路径衰减损耗===============================//
double L_duoyuan_zi_lu_jing(double ph[], int N, int u, int v, double f, double lined, double ht, double hr, int y3, int y4, double a, int e, double re)
{
    //u是子路径的起点，v是子路径的终点。
    double lamda, dd, L_dB;
    int p, q, m;
    int i;     
    m = v - u + 1;
    lamda = 300.0 / f;
    dd = lined / (N - 1);
    p = u + 1; q = v - 1;

    for (i = 1; i < m; i++)
    {
        if (p < v && ph[p - 1] > ph[p])
            p = p + 1;

        if (q > u && ph[q - 1] > ph[q - 2])
            q = q - 1;
    }

    double duv, dui, div, hri, hti, Cf, Ah;
    double *hzi;
	hzi = new double[m - 2];
    double *F;
	F = new double[m - 2];
    double *c;
	c = new double[m - 2];
    double *Li;
	Li = new double[m - 2];
   
    if (p == q)
        L_dB = 0;
    else
    {
        duv = dd * (v - u);
        if (m > 2)
        {
            for (i = 1; i < m - 1; i++)
            {
                dui = dd * i;
                div = duv - dui;
                F[i - 1] = sqrt(lamda * dui * div / duv);
                hri = (ph[u - 1] * div + ph[v - 1] * dui) / duv;
                hti = ph[u + i - 2] + dui * div / (2 * re);
                hzi[i - 1] = hri - hti;
                c[i - 1] = hzi[i - 1] / F[i - 1];
                
            }

            Cf = c[0];
            for (i = 1; i < m - 2; i++)
            {
                Cf = min(Cf, c[i]);
            }
           

            Ah = OutVisual(ph, N, lined, f, ht, hr, y3, y4, a, e, re);

            

            for (i = 0; i < m - 2; i++)
            {
                Li[i] = 10 * log10(1 - 5 / 3 * Cf / F[i]) + Ah;
               
            }
            L_dB = Li[0];
            for (i = 1; i < m - 2; i++)
            {
                L_dB = max(Li[i], L_dB);
            }
            if (L_dB < 0)

                L_dB = 0;
        }

        else
            L_dB = 0;

    }

	delete[]hzi;
	delete[]F;
	delete[]c;
	delete[]Li;

    return L_dB;
}


//===========================多个圆形障碍物的衰减计算==============================//
double duoyuan_u_L(double ph[], double phn[], int N, double u, double f, int nmax[], int m[], double lined, double ht, double hr, int y3, int y4, double a, int e, double re)
{
    double L_dB;
    int n1, n2, nmin;
    int m1, m2, i;
    n1 = n2 = m1 = m2 = 0;
    nmin = nmax[0];
    for (i = 1; i < 3; i++)
    {
        nmin = min(nmin, nmax[i]);
    }
    if (nmin == 0)
    {
        if (nmax[0] == 0)
        {
            n1 = nmax[1] + 1;
            n2 = nmax[2] + 1;
            m1 = m[1]; m2 = m[2];
        }
        else if (nmax[2] == 0)
        {
            n1 = nmax[0] + 1;
            n2 = nmax[1] + 1;
            m1 = m[0]; m2 = m[1];
        }

        double lined1, lined2;
       
        lined1 = lined / (N - 1) * (n2 - 1); lined2 = lined / (N - 1) * (N - n1);
        double *ph1;
		ph1 = new double[n2];
        double *ph2;
		ph2 = new double[N - n1 + 1];
        for (i = 0; i < n2; i++)
            ph1[i] = ph[i];
        for (i = 0; i < N - n1 + 1; i++)
            ph2[i] = ph[i + n1 - 1];
        double hn1max, hn2max, L1, L2, Rm1, Rm2;
        FUN FeiNieRm1;
        Rm1 = FeiNieRm1.fei_nie_R(ph, n2, f, lined1, n1);
        hn1max = phn[n1 - 2];
        L1 = yuan_u_L(u, Rm1, hn1max);


        Rm2 = FeiNieRm1.fei_nie_R(ph2, N - n1 + 1, f, lined2, n2 - n1 + 1);
        hn2max = phn[n2 - 2];
        L2 = yuan_u_L(u, Rm2, hn2max);

       
        ///////////////////////////////////////////////
        int x1, y1, x2, y2;
        double Lw1x1, Ly1z1, Ly2z2;
        double s11, s12, s21, s22, Pa, Pb, Cn;
        x1 = n1 - m1; y1 = n1 + m1;
        x2 = n2 - m2; y2 = n2 + m2;
        if (x1 - 1 > 1)
            Lw1x1 = L_duoyuan_zi_lu_jing(ph, N, 1, x1, f, lined, ht, hr, y3, y4, a, e, re);
        else
            Lw1x1 = 0;
        if (n2 - y1 > 1)
            Ly1z1 = L_duoyuan_zi_lu_jing(ph, N, y1, n2, f, lined, ht, hr, y3, y4, a, e, re);
        else
            Ly1z1 = 0;
        if (N - y2 > 1)
            Ly2z2 = L_duoyuan_zi_lu_jing(ph, N, y2, N, f, lined, ht, hr, y3, y4, a, e, re); //这里有问题，为什么前面两个没有问题，而只有这个有问题呢？？？？？？？？？//==这个问题已经解决了，因为问题出现在
        else                                                                  //N-y2<0,系统无法为指针找到内存位置
            Ly2z2 = 0;
        s11 = (n1 - 1) * lined / (N - 1); s12 = (n2 - n1) * lined / (N - 1);
        s21 = s12; s22 = (N - n2) * lined / (N - 1);
        Pa = s11 * (s12 * (s11 + s12 + s22)) * (s22 * (s11 + s12 + s22));
        Pb = s11 * s22 * (s11 + s12) * (s21 + s22);
        Cn = sqrt(Pa / Pb);
        
        L_dB = L1 + L2 + Lw1x1 + Ly1z1 + Ly2z2 - 20 * log10(Cn);
        if (L_dB < 0)
            L_dB = 0;
        else
            L_dB = L_dB + 0.001;

		delete[]ph1;
		delete[]ph2;

    }
    else
    {
        int n3;
        double lined1, lined2, lined3;
        n1 = nmax[0] + 1; n2 = nmax[1] + 1; n3 = nmax[2] + 1;
        double *ph1;
		ph1 = new double[n2];
        double *ph2;
		ph2 = new double[n3 - n1 + 1];
        double *ph3;
		ph3 = new double[N - n2 + 1];
        for (i = 0; i < n2; i++)
            ph1[i] = ph[i];
        for (i = 0; i < n3 - n1 + 1; i++)
            ph2[i] = ph[i + n1 - 1];
        for (i = 0; i < N - n2 + 1; i++)
            ph3[i] = ph[i + n2 - 1];

        double hn1max, hn2max, hn3max, L1, L2, L3, Rm1, Rm2, Rm3;


        lined1 = lined / (N - 1) * (n2 - 1);
        lined2 = lined / (N - 1) * (n3 - n1);
        lined3 = lined / (N - 1) * (N - n2);
        FUN FeiNieRm1;
        Rm1 = FeiNieRm1.fei_nie_R(ph1, n2, f, lined1, n1);
        hn1max = phn[n1 - 2];
        L1 = yuan_u_L(u, Rm1, hn1max);


        Rm2 = FeiNieRm1.fei_nie_R(ph2, n3 - n1 + 1, f, lined2, n2 - n1 + 1);
        hn2max = phn[n2 - 2];
        L2 = yuan_u_L(u, Rm2, hn2max);

        Rm3 = FeiNieRm1.fei_nie_R(ph3, N - n2 + 1, f, lined3, n3 - n2 + 1);
        hn3max = phn[n3 - 2];
        L3 = yuan_u_L(u, Rm3, hn3max);


        int x1, y1, x2, y2, x3, y33;
        double Lw1x1, Ly1z1, Ly2z2, Ly3z3;
        double s11, s12, s21, s22, s31, s32, Pa, Pb, Cn;

        x1 = n1 - m[0]; y1 = n1 + m[0]; x2 = n2 - m[1]; y2 = n2 + m[1]; x3 = n3 - m[2]; y33 = n3 + m[2];
        if (x1 - 1 > 1)
            Lw1x1 = L_duoyuan_zi_lu_jing(ph, N, 1, x1, f, lined, ht, hr, y3, y4, a, e, re);
        else
            Lw1x1 = 0;
        if (n2 - y1 > 1)
            Ly1z1 = L_duoyuan_zi_lu_jing(ph, N, y1, n2, f, lined, ht, hr, y3, y4, a, e, re);
        else
            Ly1z1 = 0;
        if (n3 - y2 > 1)
            Ly2z2 = L_duoyuan_zi_lu_jing(ph, N, y2, n3, f, lined, ht, hr, y3, y4, a, e, re);
        else
            Ly2z2 = 0;
        if (N - y33 > 1)
            Ly3z3 = L_duoyuan_zi_lu_jing(ph, N, y3, N, f, lined, ht, hr, y3, y4, a, e, re);
        else
            Ly3z3 = 0;

        s11 = (n1 - 1) * lined / (N - 1); s12 = (n2 - n1) * lined / (N - 1);
        s21 = s12; s22 = (n3 - n2) * lined / (N - 1);
        s31 = s22; s32 = (N - n3) * lined / (N - 1);
        Pa = s11 * (s12 * (s11 + s12 + s22 + s32)) * (s22 * (s11 + s12 + s22 + s32)) * (s32 * (s11 + s12 + s22 + s32));
        Pb = s11 * s32 * (s11 + s12) * (s21 + s22) * (s31 + s32);
        Cn = sqrt(Pa / Pb);

        L_dB = L1 + L2 + L3 + Lw1x1 + Ly1z1 + Ly2z2 + Ly3z3 - 20 * log10(Cn);
        if (L_dB < 0)
            L_dB = 0;
        else
            L_dB = L_dB + 0.001;

		delete[]ph1;
		delete[]ph2;
		delete[]ph3;

    }
    return L_dB;
}


double chazhi(double x, double y, double **Grid, int row, int column)
{
	double za, zb, zc, zd, a, b, h;

	if (x >= row - 1 || y >= column-1)
	{
			if (x >= row - 1 && y >= column - 1)
				h = Grid[(int)(floor(x))][(int)(floor(y))];
			else if (x >= row - 1)
				h = Grid[(int)(floor(x))][(int)(round(y))];
			else 
				h = Grid[(int)(round(x))][(int)(floor(y))];

	}
	else if (x <= 0 || y <= 0)
	{
		if (x <= 0 && y <= 0)
			h = Grid[(int)(ceil(x))][(int)(ceil(y))];
		else if (x <= 0)
			h = Grid[(int)(ceil(x))][(int)(round(y))];
		else
			h = Grid[(int)(round(x))][(int)(ceil(y))];

	}
	else
	{
		
		za = Grid[(int)(floor(x))][(int)(ceil(y))];
		zb = Grid[(int)(ceil(x))][(int)(ceil(y))];
		zc = Grid[(int)(ceil(x))][(int)(floor(y))];
		zd = Grid[(int)(floor(x))][(int)(floor(y))];
		a = x - (float)floor(x);
		b = y - (float)floor(y);
		h = (1 - a) * (1 - b) * zd + a * (1 - b) * zc + a * b * zb + (1 - a) * b * za;
	}
	return h;
}

double udistinguish(int N, double xt_ce, double yt_ce, double xr_ce, double yr_ce, double **Grid, double *h0, int hymax, double dd, double angle2, double lined, int row, int column)
{
	double *yy = new double[N];
	yy[0] = chazhi(xt_ce, yt_ce, Grid, row, column);
	yy[N - 1] = chazhi(xr_ce, yr_ce, Grid, row, column);

	for (int i = 1; i < N - 1; i++)
	{
		yy[i] = chazhi(xt_ce + (xr_ce - xt_ce) / (N - 1)*i, yt_ce + (yr_ce - yt_ce) / (N - 1)*i, Grid, row, column);
	}


	for (int i = 0; i < N; i++)
	{
		if (yy[i] < h0[i])
			yy[i] = 0;
	}
	int P = 0, Q = 0;
	for (int i = hymax; i > 0; i--)
	{
		if (yy[i] > 0)
			P = i;
		else
		{
			P = i - 1;
			break;
		}
	}
	for (int i = hymax + 2; i < N; i++)
	{
		if (yy[i] > 0)
			Q = i;
		else
		{
			Q = i - 1;
			break;
		}
	}
	double ro, X, u2;
	ro = dd*(Q - P) / (cos(angle2));
	X = hymax*dd / lined; //水平距离相除
	u2 = 2.02 * pow(X * (1 - X) * lined / ro, 2.0 / 3.0);
	return u2;
}


double L_danfengcemian(double acc, double lined, double f, double ht, double hr, double hymax, double hyn_max)
{
	double lamda = 300 / f;
	double dan, dnb, d1, d2, vo, Jvp;
	dan = acc*(hymax - 1);
	dnb = lined - dan;
	d1 = sqrt(hyn_max*hyn_max + dan*dan); //障碍物顶点到发射点的距离
	d2 = sqrt(hyn_max*hyn_max + dnb*dnb); //障碍物顶点到接受的距离
	vo = 0.0316*hyn_max*sqrt(2 * (d1 + d2) / (lamda*d1*d2) * 1000);
	if (vo > -0.78)
		Jvp = 6.9 + 20 * log10(sqrt((vo - 0.1) * (vo - 0.1) + 1) + vo - 0.1);
	else
		Jvp = 0;
	return Jvp;
}

double cemianzhangai(int N, double *h, double ht, double hr, double lined, double xt, double yt, double xr, double yr, double f, double dx1, double dy1, double dx, double dy, double **stt_left, double **stt_right, double *stt_lefth, double *stt_righth, double **Grid, double *ddd, double d1, double d2, double dab, double acc, double angle1, double angle2, double hh, int y3, int y4, double aa, int e, double re, int row, int column)
{
	double *hnn1 = new double[N - 2];//上余隙
	double *hnn2 = new double[N - 2];//下余隙
	double *h0 = new double[N];//天线高度
	double L_dB, L_dB1, L_dB2;
	for (int i = 0; i < N; i++)
	{
		h0[i] = h[0] + ht + (h[N - 1] + hr - h[0] - ht) * (ddd[i] / lined); //天线顶端连线 高程值
	}
	double y0, y00, y01, xt1, yt1, xr1, yr1, lined11;
	y0 = 0; y00 = 0; y01 = 0;//侧面障碍参数
	for (int i = 1; i < N-1; i++)//左侧面
	{
		if (stt_lefth[i] >= h0[i])
		{
			xt1 = xt + i * dx1;
			yt1 = yt + i * dy1;
			xr1 = stt_left[i][0];
			yr1 = stt_left[i][1];
			//两点剖面
			lined11 = sqrt(((xr1 - xt1) * dx) * ((xr1 - xt1) * dx) + ((yr1 - yt1) * dy) * ((yr1 - yt1) * dy));
			double *h11 = new double[N];
			double dx11, dy11;
			double st1[2];

			dx11 = (xr1 - xt1) / (N - 1);
			dy11 = (yr1 - yt1) / (N - 1);
			for (int ii = 0; ii < N; ii++)
			{
				st1[0] = xt1 + ii * dx11;
				st1[1] = yt1 + ii * dy11;
				h11[i] = chazhi(st1[0], st1[1], Grid, row, column);

			}
			for (int m = 0; m < N; m++)
			{
				if (h11[m] > h0[i])
				{
					hnn1[i-1] = -lined11 / (N - 1) * m;
					y00 = 1;
					break;
				}
				else
					hnn1[i-1] = -999;
			}
		}
		else
			hnn1[i-1] = -999;
	}

	for (int i = 1; i < N-1; i++)//右侧面
	{
		if (stt_righth[i] >= h0[i])
		{
			xt1 = xt + i  * dx1;
			yt1 = yt + i  * dy1;
			xr1 = stt_right[i][0];
			yr1 = stt_right[i][1];
			//两点剖面
			lined11 = sqrt(((xr1 - xt1) * dx) * ((xr1 - xt1) * dx) + ((yr1 - yt1) * dy) * ((yr1 - yt1) * dy));
			double *h22 = new double[N];
			double dx11, dy11;
			double st1[2];
			dx11 = (xr1 - xt1) / (N - 1);
			dy11 = (yr1 - yt1) / (N - 1);

			for (int ii = 0; ii < N; ii++)
			{
				st1[0] = xt1 + ii * dx11;
				st1[1] = yt1 + ii * dy11;
				h22[i] = chazhi(st1[0], st1[1], Grid, row, column);
			}

			for (int m = 0; m < N; m++)
			{
				if (h22[m] > h0[i])
				{
					hnn2[i-1] = -lined11 / (N - 1) * m;
					y01 = 1;
					break;
				}
				else
					hnn2[i-1] = -999;
			}
		}
		else
			hnn2[i-1] = -999;
	}
	if (y00 == 1 && y01 == 1)
		y0 = 3;
	if (y00 == 1 && y01 == 0)
		y0 = 1;
	if (y00 == 0 && y01 == 1)
		y0 = 2;

	if (y0>10000000000)  //存在侧面遮挡
	{
		double *hyn = new double[N - 2];
		double Rm, Hh0, Rm1, Rm2, hyn_max, hyn_max1, hyn_max2;
		double xt_ce, yt_ce, xr_ce, yr_ce, u2;

		if (y0 == 1)  //左侧遮挡
		{
			for (int i = 0; i < N - 2; i++)
			{
				hyn[i] = hnn1[i];
			}
			hyn_max = hyn[0];
			int hymax = 0;
			for (int i = 1; i < N - 2; i++)
			{
				hyn_max = max(hyn_max, hyn[i]);
			}
			for (int i = 0; i < N - 2; i++)
			{
				if (hyn[i] == hyn_max)
				{
					hymax = i;
					break;
				}
			}
			d1 = (hymax + 1)*acc / 1000 / cos(angle2);
			d2 = dab / 1000 - d1;
			Rm = 550 * sqrt(d1*d2 / dab / f); //nmax处的菲涅尔半径
			Hh0 = Rm / sqrt(3);
			hh = Hh0 + abs(hyn_max);
			xt_ce = xt - hh*sin(angle1) / dx;
			yt_ce = yt - hh*cos(angle1) / dy;
			xr_ce = xr - hh*sin(angle1) / dx;
			yr_ce = yr - hh*cos(angle1) / dy;
			u2 = udistinguish(N, xt_ce, yt_ce, xr_ce, yr_ce, Grid, h0, hymax, acc, angle2, lined,row,column);

			if (u2 > 3)
				L_dB = L_danfengcemian(acc, lined, f, ht, hr, hymax, hyn_max);
			else
				L_dB = yuan_u_L(u2, Rm, hyn_max);
		}
		if (y0 == 2)  //右侧遮挡
		{
			for (int i = 0; i < N - 2; i++)
			{
				hyn[i] = hnn2[i];
			}
			double hyn_max = hyn[0];
			int hymax = 0;
			for (int i = 1; i < N - 2; i++)
			{
				hyn_max = max(hyn_max, hyn[i]);
			}
			for (int i = 0; i < N - 2; i++)
			{
				if (hyn[i] == hyn_max)
				{
					hymax = i;
					break;
				}
			}
			d1 = (hymax + 1)*acc / 1000 / cos(angle2);
			d2 = dab / 1000 - d1;
			Rm = 550 * sqrt(d1*d2 / dab / f); //nmax处的菲涅尔半径
			Hh0 = Rm / sqrt(3);
			hh = Hh0 + abs(hyn_max);
			xt_ce = xt + hh*sin(angle1) / dx;
			yt_ce = yt + hh*cos(angle1) / dy;
			xr_ce = xr + hh*sin(angle1) / dx;
			yr_ce = yr + hh*cos(angle1) / dy;
			u2 = udistinguish(N, xt_ce, yt_ce, xr_ce, yr_ce, Grid, h0, hymax, acc, angle2, lined, row, column);

			if (u2 > 3)
				L_dB = L_danfengcemian(acc, lined, f, ht, hr, hymax, hyn_max);
			else
				L_dB = yuan_u_L(u2, Rm, hyn_max);
		}
		if (y0 == 3)  //两侧遮挡
		{
			int u21, u22;
			for (int i = 0; i < N - 2; i++)
			{
				hyn[i] = hnn1[i];
			}
			double hyn_max1 = hyn[0];
			int hymax = 0;
			for (int i = 1; i < N - 2; i++)
			{
				hyn_max1 = max(hyn_max1, hyn[i]);
			}
			for (int i = 0; i < N - 2; i++)
			{
				if (hyn[i] == hyn_max1)
				{
					hymax = i;
					break;
				}
			}
			d1 = (hymax + 1)*acc / 1000 / cos(angle2);
			d2 = dab / 1000 - d1;
			Rm1 = 550 * sqrt(d1*d2 / dab / f); //nmax处的菲涅尔半径
			Hh0 = Rm1 / sqrt(3);
			hh = Hh0 + abs(hyn_max1);
			xt_ce = xt - hh*sin(angle1) / dx;
			yt_ce = yt - hh*cos(angle1) / dy;
			xr_ce = xr - hh*sin(angle1) / dx;
			yr_ce = yr - hh*cos(angle1) / dy;
			u21 = udistinguish(N, xt_ce, yt_ce, xr_ce, yr_ce, Grid, h0, hymax, acc, angle2, lined, row, column);
			if (u21 > 3)
				L_dB1 = L_danfengcemian(acc, lined, f, ht, hr, hymax, hyn_max);
			else
				L_dB1 = yuan_u_L(u21, Rm1, hyn_max1);


			for (int i = 0; i < N - 2; i++)
			{
				hyn[i] = hnn2[i];
			}
			double hyn_max2 = hyn[0];
			hymax = 0;
			for (int i = 1; i < N - 2; i++)
			{
				hyn_max2 = max(hyn_max2, hyn[i]);
			}
			for (int i = 0; i < N - 2; i++)
			{
				if (hyn[i] == hyn_max2)
				{
					hymax = i;
					break;
				}
			}
			d1 = (hymax + 1)*acc / 1000 / cos(angle2);
			d2 = dab / 1000 - d1;
			Rm2 = 550 * sqrt(d1*d2 / dab / f); //nmax处的菲涅尔半径
			Hh0 = Rm2 / sqrt(3);
			hh = Hh0 + abs(hyn_max2);
			xt_ce = xt + hh*sin(angle1) / dx;
			yt_ce = yt + hh*cos(angle1) / dy;
			xr_ce = xr + hh*sin(angle1) / dx;
			yr_ce = yr + hh*cos(angle1) / dy;
			u22 = udistinguish(N, xt_ce, yt_ce, xr_ce, yr_ce, Grid, h0, hymax, acc, angle2, lined, row, column);


			if (u22 > 3)
				L_dB2 = L_danfengcemian(acc, lined, f, ht, hr, hymax, hyn_max);
			else
				L_dB2 = yuan_u_L(u22, Rm2, hyn_max2);

			L_dB = L_dB1 + L_dB2;
		}

	}
	else  //不存在侧面遮挡
	{
		double h1, h2, hm, dlos;
		h1 = h[0] + ht;
		h2 = h[N - 1] + hr;
		hm = sqrt(h1) + sqrt(h2);
		dlos = sqrt(2.0 * re) * hm;
		if (lined >= dlos)
			L_dB = OutVisual(h, N, lined, f, ht, hr, y3, y4, aa, e, re);

		else
			L_dB = InVisual(h, N, lined, f, ht, hr, y3, y4, aa, e, re);
	}

	return L_dB;
}


double ITUrp::ITUraoshe(double xt, double yt, double xr, double yr, double **Grid, int row, int column, double acc, double f, double ht, double hr, int y3, int y4, double aa, int e, double re, double dx, double dy)
{
	int N;
	double L_dB = 0.0, L_dB1 = 0.0, L_dB2 = 0.0;
	double lined;
	lined = sqrt(((xr - xt) * dx) * ((xr - xt) * dx) + ((yr - yt) * dy) * ((yr - yt) * dy));
	N = (int)(lined / acc) + 1;
	double *h = new double[N];
	double dx1, dy1;
	double st[2];

	dx1 = (xr - xt) / (N - 1);
	dy1 = (yr - yt) / (N - 1);

	int i;
	for (i = 0; i < N; i++)
	{
		st[0] = xt + i * dx1;
		st[1] = yt + i * dy1;
		h[i] = chazhi(st[0], st[1], Grid, row, column);
	}

	/////////////////////////////////////////////////////////////////////////////
	double dan, dnb;
	double* hn;
	hn = new double[N - 2];
	for (i = 1; i < N - 1; i++)
	{
		dan = lined / (N - 1) * i;  //到发射点的距离
		dnb = lined - dan;      //到接受点的距离
		hn[i - 1] = h[i] + (dan * dnb / (2 * re)) - (((h[0] + ht) * dnb + (h[N - 1] + hr) * dan) / lined);  //菲涅尔余隙，不是菲涅尔半径

	}

	double y1;
	double hya, hyb, dab, Rmax, hh;

	hya = h[0] + ht; hyb = h[N - 1] + hr;
	dab = sqrt(lined * lined + (hyb - hya) * (hyb - hya));
	Rmax = 550 * sqrt((dab / f / 4000)); //第一菲涅尔半径最大值,d1=d2=dab/2时


	double hmax, hmin;
	hmax = h[0]; hmin = h[0];
	for (i = 1; i < N; i++)
	{

		hmax = max(hmax, h[i]);
		hmin = min(hmin, h[i]);
	}
	hh = 0.8 * (hmax - hmin); //地形不规则度

	if (hh <= 0.1 * Rmax) //地球是否光滑的判断条件
		y1 = 1;        //y1=1表示地形为光滑的
	else
		y1 = 0;
	//y1=1代表地球为光滑的，否则可能有障碍物
	if (y1 == 1)
	{
		double* ddd = new double[N];
		for (i = 0; i < N; i++)
		{
			ddd[i] = acc*i;   //发射接收点之间间隔取点
		}
		double d1, d2;
		double angle1 = atan(abs((yr - yt) / (xr - xt)));
		double angle2 = atan(abs((hyb - hya) / lined));
		double* R;
		R = new double[N - 2];
		for (i = 1; i < N - 1; i++)
		{
			d1 = i*acc / 1000 / cos(angle2);   //天线顶端连线间隔取点
			d2 = dab / 1000 - d1;
			R[i-1] = 550 * sqrt(d1*d2 / (dab /1000) / f);
		}

		double **stt_left;
		stt_left = new double*[N];
		for (int i = 0; i < N; i++)
		{
			stt_left[i] = new double[2];
		}

		double **stt_right;
		stt_right = new double*[N];
		for (int i = 0; i < N; i++)
		{
			stt_right[i] = new double[2];
		}

		stt_left[0][0] = xt;
		stt_left[0][1] = yt;
		stt_left[N - 1][0] = xr;
		stt_left[N - 1][1] = yr;
		double* stt_lefth;
		stt_lefth = new double[N];
		stt_lefth[0] = hya;
		stt_lefth[N-1] = hyb;

		stt_right[0][0] = xt;
		stt_right[0][1] = yt;
		stt_right[N - 1][0] = xr;
		stt_right[N - 1][1] = yr;
		double* stt_righth;
		stt_righth = new double[N];
		stt_righth[0] = hya;
		stt_righth[N-1] = hyb;
		for (i = 1; i < N - 1; i++)
		{
			d1 = i*acc / 1000 / cos(angle2);
			d2 = dab / 1000 - d1;
			stt_left[i][0] = xt + (xr - xt) / (N - 1)*i - R[i-1] / sqrt(3)*sin(angle1) / dx;
			stt_left[i][1] = yt + (yr - yt) / (N - 1)*i - R[i-1] / sqrt(3)*cos(angle1) / dy;
			stt_right[i][0] = xt + (xr - xt) / (N - 1)*i + R[i-1] / sqrt(3)*sin(angle1) / dx;
			stt_right[i][1] = yt + (yr - yt) / (N - 1)*i + R[i-1] / sqrt(3)*cos(angle1) / dy;
			stt_lefth[i] = chazhi(stt_left[i][0], stt_left[i][1], Grid, row, column);
			stt_righth[i] = chazhi(stt_right[i][0], stt_right[i][1], Grid, row, column);
		}

		L_dB = cemianzhangai(N, h, ht, hr, lined, xt, yt, xr, yr, f, dx1, dy1, dx, dy, stt_left, stt_right, stt_lefth, stt_righth, Grid, ddd, d1, d2, dab, acc, angle1, angle2, hh, y3, y4, aa, e, re, row, column);
		delete []ddd;
		delete []R;
		for (int i = 0; i < N; i++)
		{
			delete[] stt_left[i];
		}
		delete[] stt_left;
		for (int i = 0; i < N; i++)
		{
			delete[] stt_right[i];
		}
		delete[] stt_right;
		delete[]stt_lefth;
		delete[]stt_righth;
	}
    else
    {
        int ObstacleNumber = 0;
        double d[2];
        int pointnum;

        //d[0] = fabs(xt - xr);
        //d[1] = fabs(yt - yr);

      //  pointnum = (int)(Math.Ceiling((Math.Sqrt(d[0] * d[0] + d[1] * d[1])) / 10) * 10 - 1);
        pointnum = 2 * N;  //判定障碍物时精度需要更高？
        double *y;
		y = new double[pointnum + 2];
        
        for (i = 0; i < pointnum + 2; i++)
        {
			st[0] = xt + i * (xr - xt) / (pointnum + 1);
			st[1] = yt + i * (yr - yt) / (pointnum + 1);
			y[i] = chazhi(st[0], st[1], Grid, row, column);


        }
       
        ////////////////////////////////////////////////////////////////////////////////
        Obstacle UDistinguish;
        UDistinguish.U_Distinguish(y, lined, pointnum, f, ht, hr, re);//这里如果采样点越密精确度越高，U值越精确，但是计算量会增大。
     
        Obstacle OBstNUM;
        OBstNUM.ObstacleNum(hn, h, lined, N, f);

        ObstacleNumber = OBstNUM.y5;

		if (ObstacleNumber == 0)  //障碍物个数为0，不存在孤立障碍物
		{
			L_dB = 0;
		}
        else //存在障碍物
        {
            if (UDistinguish.u > 3) //障碍物类型为刀刃
            {
				if (ObstacleNumber == 1)
				{
					L_dB = L_danfeng(h, hn, N, lined, f, re);//障碍物为单个刀刃障碍物
				}
                if (ObstacleNumber == 2)
                {
                    int num[3];    //三个障碍物中心位置参数数组
                    num[0] = OBstNUM.Nmax[0];//第一个障碍物中心位置
                    num[1] = OBstNUM.Nmax[1];//第二个障碍物中心位置
                    num[2] = OBstNUM.Nmax[2];//第三个障碍物中心位置
                    L_dB = L_doublefeng(h, hn, N, lined, f, re, num);//障碍物为双重刀刃障碍物
                }
                if (ObstacleNumber == 3)
                {
                    L_dB = L_duofeng(h, hn, N, lined, f, ht, hr, re);//障碍物为三个刀刃障碍物
                }
            }
            else  //障碍物类型为圆形
            {
                
                L_dB = yuan_u_L(UDistinguish.u, UDistinguish.Rm, UDistinguish.ynmax);
               
            }
        }
		delete []y;

    }

	delete []h;
	delete []hn;
	
    L_dB = L_dB + 32.45 + 20 * log10(lined / 1000) + 26 * log10(f);
    return L_dB;
}
//======================================================================================================//
