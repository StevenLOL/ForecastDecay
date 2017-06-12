// attenuation_ra.cpp: implementation of the attenuation_ra class.
//
//////////////////////////////////////////////////////////////////////

#include "attenuation_ra.h"
#include "math.h"
#include <stdio.h>
#define PI 3.1415927

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

attenuation_ra::attenuation_ra()
{

}

attenuation_ra::~attenuation_ra()
{

}
static double lamde;		// 微波的波长（km）
static double h1;			// 发射天线的高度
static double h2;			// 计算出来的接收天线的高度
static double Hts;			// 发射端天线的平均高度
static double Hrs;			// 接收端天线的平均高度
static double maxval;		// 发射端的物理视界仰角
static int maxnum;			// 发射端视界临界点在路径上的位置
static double maxval1;		// 接收端的物理视界仰角
static int maxnum1;			// 接收端视界临界点在路径上的位置
static int type;			// 视界or非视界的标志位

static double thetat, thetar;	// 发射端（接收端）视界俯仰角
static double theta;			// 路径的角向距离
static int ilt, ilr;			// 发射端（接收端）视界临界点在链路上的位置
static double dlt , dlr;		// 从发射天线和接收天线到相应的视界的距离。
static double dct , dcr;		// 从发射端（接收端）到沿大圆路径上跨越陆地的距离

static double Hte;				// 发射端天线有效高度
static double Hre;				// 接收端天线的有效高度
static double hm ;				// 地面粗糙度
// = = = = = = = = = = = = = = = 视距 = = = = = = = = = = = = = = = = = = = = = = = = //  
static double L_b0p;			// 视距传播P%时间内不超过的基本损耗
static double L_b0b;			// 视距传播b0%时间内不超过的基本损耗
static double Ag;				// 大气吸收损耗

// = = = = = = = = = = = = = = = 绕射 = = = = = = = = = = = = = = = = = = = = = = = = //
static double Ld50;				// 绕射损耗中值Ld50√
static double Ldb;				// b%时间内不超过的绕射损耗Ldβ √
static double L_dp;				// 对于p%时间不超过的绕射损耗√
static  double L_bd50;			// 与绕射有关的基本传输损耗中值√
static  double L_bd;			// 与绕射有关的对于p%时间不超过的基本传输损耗√
static double Fi;				// √
void attenuation_ra::shiju(double f,double p,double b0, double d, double dlt, double dlr, double Ag)
{
	double L_bfsg; // 视距传播的基本传输损耗
	double Esp,Esb;
	L_bfsg = 92.5 + 20 * log10(f) + 20 *log10(d) +Ag;
	Esp = 2.6 * (1.0 - exp(-0.1*(dlt+dlr))) * log10(p/50.0); // 在p%时间内由多径和聚焦效应引起的校正项
	Esb = 2.6 * (1.0 - exp(-0.1*(dlt+dlr))) * log10(b0/50.0); // 在b0%时间内由多径和聚焦效应引起的校正项
	L_b0p = L_bfsg + Esp; // 在p%时间内不超过的由视距传播引起的基本传播损耗
	L_b0b = L_bfsg + Esb; // 在b0%时间内不超过的由视距传播引起的基本传播损耗
}

double calcuI(double x)
{
	// 反累积正太函数：基于绕射损耗成对数正太分布的内插因子
    double I;
    double sigema;
    double T;
    double C0 = 2.515516698, C1 = 0.802853, C2 = 0.010328;
    double D1 = 1.432788, D2 = 0.189269, D3 = 0.001308;

    T = sqrt(-2.0 * log(x));
    sigema = (((C2 * T + C1) * T) + C0) / (((D3 * T + D2) * T + D1) * T + 1);
    I = sigema - T;

    return I;
}

double calcuJ(double v)
{
	// 单楔形绕射损耗近似值的函数
    double J;
    if (v > -0.78)
        J = 6.9 + 20 * log10(sqrt((v - 0.1) * (v - 0.1) + 1) + v - 0.1);
    else J = 0;
    return J;
}

void attenuation_ra::raosheParam(int pointnum, double di[], double hi[], double d, double ae, double abeite)
{
    double keseim = 0;  //用于所有路径斜率的校正项ζm
    double* Hi = new double[pointnum + 1];  //垂直净空Hi
    double* vmi50 = new double[pointnum + 1];
    double vm50 = -1000;    //绕射参数
    int im50 = 0;   //具有最大值vm50的剖面点的指数

    keseim = cos(atan(1e-3 * (Hrs - Hts) / d));
    // ============================ 计算曲面路径与天线射线间的最小间距============================ //
    for (int i = 0; i < pointnum + 1; i++)
    {
        Hi[i] = hi[i] + 1e3 * di[i] * (d - di[i]) / (2 * ae) - (Hts * (d - di[i]) + Hrs * di[i]) / d;
        vmi50[i] = keseim * Hi[i] * sqrt(2e-3 * d / (lamde * di[i] * (d - di[i])));

        if (vm50 < vmi50[i])
        {
            vm50 = vmi50[i];
            im50 = i;
        }
    }
	// ================================ 计算绕射损耗中值 ========================================== //
    double Lm50 = 0;    // 主要边界的楔形绕射损耗中值Lm50
    double Lt50 = 0;    // 发射机测次边界的楔形绕射损耗中值Lt50
    double Lr50 = 0;    // 接收机侧次边界的楔形绕射损耗中值Lr50

    double keseit = cos(atan(1e-3 * (Hi[im50] - Hts) / di[im50]));          // 校正项
    double keseir = cos(atan(1e-3 * (Hrs - Hi[im50]) / (d - di[im50])));    // 校正项
    int it50 = 0;   // 发射机侧次边界剖面点的指数
    int ir50 = 0;   // 接收机侧次边界剖面点的指数

    Lm50 = calcuJ(vm50);
    // 若Lm50 = 0，则绕射损耗中值Ld50和Ldβ均为0，不必再计算绕射。
    if (Lm50 == 0)
    {
        Ld50 = 0;
        Ldb = 0;
    }
    // 若Lm50 ≠ 0，则如下计算主边界的发射机和接收机侧的次边界可能引起的绕射损耗
    else
    {
        //------------------------- 发射机测次边界的绕射损耗中值Lt50 ----------------------------- //
        if (im50 == 1)   // 不存在发射机侧次边界，相关绕射损耗Lt50置零
            Lt50 = 0;
        else
        {
            double* vti50 = new double[pointnum + 1];
            double vt50 = -1000;    // 绕射参数

            for (int i = 1; i < im50; i++)
            {
                Hi[i] = hi[i] + 1e3 * di[i] * (di[im50] - di[i]) / (2 * ae) - (Hts * (di[im50] - di[i]) + hi[im50] * di[i]) / di[im50];
                vti50[i] = keseit * Hi[i] * sqrt(2e-3 * di[im50] / (lamde * di[i] * (di[im50] - di[i])));
                if (vt50 < vti50[i])
                {
                    vt50 = vti50[i];
                    it50 = i;
                }
            }
			delete []vti50;
            Lt50 = calcuJ(vt50);
        }
        //----------------------------- 接收机侧次边界的绕射损耗中值Lr50 ----------------------------- //
        if (im50 == pointnum - 1)
            Lr50 = 0;
        else
        {
            double* vri50 = new double[pointnum + 1];
            double vr50 = -1000;    //绕射参数

            for (int i = im50 + 1; i < pointnum; i++)
            {
                Hi[i] = hi[i] + 1e3 * (di[i] - di[im50]) * (d - di[i]) / (2 * ae) - (hi[im50] * (d - di[i]) + Hrs * (di[i] - di[im50])) / (d - di[im50]);
                vri50[i] = keseir * Hi[i] * sqrt(2e-3 * (d - di[im50]) / (lamde * (di[i] - di[im50]) * (d - di[i])));
                if (vr50 < vri50[i])
                {
                    vr50 = vri50[i];
                    ir50 = i;
                }
            }
            Lr50 = calcuJ(vr50);
			delete []vri50;
        }
        Ld50 = Lm50 + (1 - exp(-Lm50 / 6)) * (Lt50 + Lr50 + 10 + 0.04 * d);

        // 4.2.2 在β0% 时间内不超过的绕射损耗Ldβ
        // 采用有效地球半径 aβ来计算
        // 在β0%时间内不超过的主边界楔形绕射损耗Lmβ
        double Lmbita;
        double vmbita;  // 绕射参数vmβ
        double Himbita = hi[im50] + 1e3 * di[im50] * (d - di[im50]) / (2 * abeite) - (Hts * (d - di[im50]) + Hrs * di[im50]) / d;

        vmbita = keseim * Himbita * sqrt(2e-3 * d / (lamde * di[im50] * (d - di[im50])));
        Lmbita = calcuJ(vmbita);

        // 在β0%时间内不超过的发射机侧次边界楔形绕射损耗Ltβ
        double Ltbita;

        // 若 Lt50 = 0，Ltβ将为0；反之，计算Ltβ如下：
        if (Lt50 == 0)
            Ltbita = 0;
        else
        {
            double vtbita;
            double Hitbita;

            Hitbita = hi[it50] + 1e3 * di[it50] * (di[im50] - di[it50]) / (2 * abeite) - ( Hts * (di[im50] - di[it50]) + hi[im50] * di[it50]) / di[im50];
            vtbita = keseit * Hitbita * sqrt(2e-3 * di[im50] / ( lamde * di[it50] * (di[im50] - di[it50])));
            Ltbita = calcuJ(vtbita);
        }

        // 在β0%时间内不超过的接收机测次边界楔形绕射损耗Lrβ
        double Lrbita;

        if (Lr50 == 0)
            Lrbita = 0;
        else
        {
            double vrbita;
            double Hirbita;

            Hirbita = hi[ir50] + 1e3 * (di[ir50] - di[im50]) * (d - di[ir50]) / (2 * abeite) - (hi[im50] * (d - di[ir50]) + Hrs * (di[ir50] - di[im50])) / (d - di[im50]);
            vrbita = keseir * Hirbita * sqrt(2e-3 * (d - di[im50]) / ( lamde * (di[ir50] - di[im50]) * (d - di[ir50])));
            Lrbita = calcuJ(vrbita);
        }
        // 在β0%时间内不超过的各种边界损耗组合在一起：
        Ldb = Lmbita + (1 - exp(-Lmbita / 6)) * (Ltbita + Lrbita + 10 + 0.04 * d);
    }
	
	
	delete []vmi50;
	delete []Hi;
}

void attenuation_ra::raoshe(double f, double p, double beite0, double omiga, double d, int pointnum, double di[], double hi[], double ae, double abeite, double Ag)
{
    // 由绕射模型计算总损耗需要下列3个值：
    // 对于p%时间不超过的绕射损耗 L_dp
    // 与绕射有关的基本传输损耗中值 L_bd50
    // 与绕射有关的对于p%时间不超过的基本传输损耗 L_bd 

    double L_bfsg = 0;
    raosheParam(pointnum, di, hi, d, ae, abeite);

    L_bfsg = 92.5 + 20 * log10(f) + 20 * log10(d) + Ag;

    if (p == 50)
        Fi = 0;
    else if (p < 50 && p > beite0)
        Fi = calcuI(p / 100) / calcuI(beite0 / 100);
    else
        Fi = 1;

    L_dp = Ld50 + Fi * (Ldb - Ld50); // 对于p%时间不超过的绕射损耗 L_dp
    L_bd50 = L_bfsg + Ld50;          // 与绕射有关的基本传输损耗中值 L_bd50
    L_bd = L_b0p + L_dp;             // 与绕射有关的对于p%时间不超过的基本传输损耗 L_bd 
}

double calcufai(double rp, double rt, double a, double b, double c, double d)
{
    double result;

    result = pow(rp, a) * pow(rt, b) * exp(c * (1 - rp) + d * (1 - rt));

    return result;
}

double calcuexpr(double f, double yita, double rt, double a, double b, double c, double d)
{
    double result;

    result = a * yita * exp(b * (1 - rt)) / pow((f - c), 2) + d * yita * yita;

    return result;
}

double calcug(double f, double fi)
{
    double g;

    g = 1 + pow((f - fi) / (f + fi), 2);

    return g;
}

double AtmoGases(double f, double d, double omiga)
{
    // 从海平面到10 km高度范围内，由干燥空气和水汽造成的无线电波特征衰减
    double Ag = 0;  // 总的气体吸收，返回值
    double ro = 0;  // 大气衰减γo(dB/km)
    double rw = 0;  // 对在水汽中的衰减，大气衰减值γw (dB/km)

    // 中间参数
    double rt = 0, rp = 0;
    double kesei1 = 0, kesei2 = 0, kesei3 = 0;  // 计算 ro 中间参数
    double yita1 = 0, yita2 = 0;  // 计算 rw 中间参数

    double pres = 1013.25; // 气压（hPa）
    double temp = 15;   // 温度（°C）
    double rou = 7.5 + 2.5 * omiga;    // 水蒸气密度(g/m^3)

    rt = 288 / (273 + temp);
    rp = pres / 1013;

    kesei1 = calcufai(rp, rt, 0.0717, -1.8132, 0.0156, -1.6515);
    kesei2 = calcufai(rp, rt, 0.5146, -4.6368, -0.1921, -5.7416);
    kesei3 = calcufai(rp, rt, 0.3414, -6.5851, 0.2130, -8.5854);

    yita1 = 0.955 * rp * pow(rt, 0.68) + 0.006 * rou;
    yita2 = 0.735 * rp * pow(rt, 0.5) + 0.0353 * pow(rt, 4) * rou;

    ro = (7.2 * pow(rt, 2.8) / (f * f + 0.34 * rp * rp * pow(rt, 1.6)) + 0.62 * kesei3 / (pow((54 - f), 1.16 * kesei1) + 0.83 * kesei2)) * f * f * rp * rp * 1e-3;
    rw = (calcuexpr(f, yita1, rt, 3.98, 2.23, 22.235, 9.42) * calcug(f, 22) + calcuexpr(f, yita1, rt, 11.96, 0.7, 183.31, 11.14) + 
		calcuexpr(f, yita1, rt, 0.081, 6.44, 321.226, 6.29) + calcuexpr(f, yita1, rt, 3.66, 1.6, 325.153, 9.22) + 
		calcuexpr(f, yita1, rt, 25.37, 1.09, 380, 0) + calcuexpr(f, yita1, rt, 17.4, 1.46, 448, 0) + 
		calcuexpr(f, yita1, rt, 844.6, 0.17, 557, 0) * calcug(f, 557) + calcuexpr(f, yita1, rt, 290, 0.41, 752, 0) * calcug(f, 752) + 
		calcuexpr(f, yita2, rt, 8.3328 * 1e4, 0.99, 1780, 0) * calcug(f, 1780)) * f * f * pow(rt, 2.5) * rou * 1e-4;
    Ag = (ro + rw) * d;
    return Ag;
}

double sanshe(double f, double p, double d, double thita, double N0, double Gt, double Gr)
{
    double L_bs;
    double Lf, Lc, Ag = 0;	// Ag为气体吸收

    Ag = AtmoGases(f, d, -1.8);
    Lf = 25 * log10(f) - 2.5 * log10(f / 2.0) * log10(f / 2.0);	// 与频率有关的损耗
    Lc = 0.051 * exp(0.055 * (Gt + Gr));
    L_bs = 190 + Lf + 20 * log10(d) + 0.573 * thita - 0.15 * N0 + Lc + Ag - 10.1 * pow((-log10(p / 50)), 0.7);
    return L_bs;
}

double attenuation_ra::atmoduct(double f, double omiga, double p, double b0, double ae, double d, double dlt, double dlr, double thita, double thitat, double thitar, double Ag)
{
    double L_ba;	// 大气波导/层反射损耗
    double A_f;     // 大气内在天线和异常传播结构之间的固定耦合损耗的综合（本地散射损耗除外）
    double A_dp;	// 在异常传播机理内与时间百分比和角度--距离有关的损耗

    // ===============================================1. 计算A_f================================================== //
    // 中间参数
    double Ast = 0, Asr = 0;	// 干扰站、被干扰站的位置屏蔽损耗
    double Act = 0, Acr = 0;	// 干扰站、被干扰站的跨海表面大气波导耦合校正量
    double thitatpp, thitarpp;

    thitatpp = thitat - 0.1 * dlt;    // mrad
    thitarpp = thitar - 0.1 * dlr;    // mrad

    // ------------------------------------------ 计算Ast，Asr ----------------------------------------- //
    if (thitatpp <= 0)
        Ast = 0;
    else
        Ast = 20 * log10(1 + 0.361 * thitatpp * sqrt(f * dlt)) + 0.264 * thitatpp * pow(f, 1 / 3);

    if (thitarpp <= 0)
        Asr = 0;
    else
        Asr = 20 * log10(1 + 0.361 * thitarpp * sqrt(f * dlr)) + 0.264 * thitarpp * pow(f, 1 / 3);

    // ------------------------------------------ 计算Act，Acr ------------------------------------------ //
    if (omiga >= 0.75 && dct <= dlt && dct <= 5)
        Act = -3 * exp(-0.25 * dct * dct) * (1 + tanh(0.07 * (50 - Hts)));
    else
        Act = 0;

    if (omiga >= 0.75 && dcr <= dlr && dcr <= 5)
        Acr = -3 * exp(-0.25 * dcr * dcr) * (1 + tanh(0.07 * (50 - Hrs)));
    else
        Acr = 0;

    A_f = 102.45 + 20 * log10(f) + 20 * log10(dlt + dlr) + Ast + Asr + Act + Acr; //!!!! dlt,dlr 

    //============================================2. 计算A_dp===============================================//
    // ------------------------------------------ 2.1 计算rd ------------------------------------------ //
    double rd;	// 比衰减
    rd = 5e-5 * ae * pow(f, 1 / 3);

    // ------------------------------------------ 2.2 计算 thitap -------------------------------------- //
    double thitap;
    double thitatp, thitarp;

    if (thitat <= 0.1 * dlt)
        thitatp = thitat;
    else
        thitatp = 0.1 * dlt;

    if (thitar <= 0.1 * dlr)
        thitarp = thitar;
    else
        thitarp = 0.1 * dlr;

    thitap = 1e3 * d / ae + thitatp + thitarp;

    // ------------------------------------------ 2.3 计算Ap ------------------------------------------ //
    double Ap;	//时间百分比p的可变性（累积分布）
    double TT;  // TT : Γ
    double bb;
    double u2, u3;
    double aa;
    double dd = 3.5;   // dd:ε
    double dlm = 0;    // 大圆路径上最长的连续陆地。正常全是陆地的情况下dlm = d
    double tal = floor(1 - exp(-4.12e-4 * pow(dlm, 2.41)));    //公式 3-133a，是否是取整函数？
    double dI = (d - dlt - dlr < 40) ? d - dlt - dlr : 40;
    double Temp = 0.0;

    aa = -0.6 - dd * 1e-9 * pow(d, 3.1) * tal;
    Temp = 500 * d * d / (ae * pow(sqrt( Hte) + sqrt( Hre), 2));
    u2 = pow(500 * d * d / (ae * pow(sqrt( Hte) + sqrt( Hre), 2)), aa);   //!!! 此系数需要调整

    u2 = (u2 > 1) ? 1 : u2;

    if ( hm <= 10)
        u3 = 1;
    else
        u3 = exp(-4.6e-5 * ( hm - 10) * (43 + 6 * dI));

    bb = b0 * u2 * u3;

    TT = 1.076 / pow(2.0058 - log10(bb), 1.012) * exp(-(9.51 - 4.8 * log10(bb) + 0.198 * log10(bb) * log10(bb)) * 1e-6 * pow(d, 1.13));
    Ap = -12 + (1.2 + 3.7e-3 * d) * log10(p / bb) + 12 * pow(p / bb, TT);
	// 异常传播机理内与时间百分数和角度—距离有关的损耗
    A_dp = rd * thitap + Ap;
	// 异常传播期间的基本传输损耗
    L_ba = A_f + A_dp + Ag;

    return L_ba;
}
  
void attenuation_ra::paramcalcu(double f, double ae, double d, double di[], double hi[], double h_radio, int pointnum, double acc,double ha, double m, double omiga)
{
    lamde = 0.3 / f;
	int i = 0;
    int n = -1;
    int hight = (int)floor(h_radio);
    h1 = 10;	// 发射天线高度
    h2 = 0;		// 接收天线高度
    double* Hci = new double[pointnum + 1];    //插值点余隙
    double Hcmin = -100;   //最小余隙
   // ==================================== 天线高度计算 ========================================== //
	// 计算传输余隙
    for (i = 1; i < pointnum; i++)
    {
        Hci[i] = ((hi[0] + hight) * (d - di[i]) + (hi[pointnum] + hight) * di[i]) / d - hi[i];
    }
	 double minHci = Hci[1];	// 最下的路径余隙
    int minHcNum = 1;			// 最小余隙点的位置
    double minHr = 0;			// 接收天线的最小高度
    for (i = 1; i < pointnum - 1; i++)
    {
        if (Hci[i] < minHci)
        {

            minHci = Hci[i];
            minHcNum = i;
        }
    }
	double tempF0 = 0.577 * sqrt(lamde * di[minHcNum] * 1000 * (d - di[minHcNum]) / d);
	while (h1 < hight)
    {
        minHr = (tempF0 - h1 - hi[0] + (d - di[minHcNum]) * di[minHcNum] / (2 * ae)+ hi[minHcNum]) * d / di[minHcNum] + h1 + hi[0] - hi[pointnum];
        if (minHr < 0)
        {
            minHr = 1;
        }
        n = 0; h2 = 0;
        if (minHr < hight)
        {
            while (h2 < hight)
            {
                h2 = (2 * n + 1) * lamde * d * 1000 / (4 * (hi[0] + h1));
                if ((h2 < hight) && (h2 > minHr))
                {
                    Hcmin = ((hi[0] + h1) * (d - di[minHcNum]) + (hi[pointnum] + h2) * di[minHcNum]) / d - hi[minHcNum];
                    if (Hcmin > -1.0 * tempF0)
                    {
                        break;
                    }
                }
                n++;
            }
        }
        if (Hcmin > -1.0 * tempF0)
        {
            break;
        }
        h1 += 1;
    }
	if (Hcmin < -1.0 * tempF0)
    {
        h1 = hight;
        h2 = hight;
    }

	// ==================================== 天线高度计算结束 ========================================== //
    Hts = hi[0] + h1;
    Hrs = hi[pointnum] + h2;

    double thitatd = (Hrs - Hts) / d - 1000 *d / (2 * ae);
    double* thitai = new double[pointnum + 1];
    double* thitaj = new double[pointnum + 1];

    thitai[0] = 0;
    thitaj[0] = (Hts - Hrs) / d - 1000 * d / (2 * ae);

    thitai[pointnum] = thitatd;
    thitaj[pointnum] = 0;

    for (i = 1; i < pointnum; i++)
    {
        thitai[i] = (hi[i] - Hts) / di[i] - 1000 * di[i] / (2 * ae);	//thitai 单位为 mrad
        thitaj[i] = (hi[i] - Hrs) / (d - di[i]) - 1000 * (d - di[i]) / (2 * ae);

        if (maxval < thitai[i])
        {
            maxval = thitai[i];
            maxnum = i;
        }
        if (maxval1 < thitaj[i])
        {
            maxval1 = thitaj[i];
            maxnum1 = i;
        }
    }

   
	if(maxval > thitatd)
	{
		type = 1;
	}
	else
	{
		type = 0;
	}
	
	if(type > 0)
    {
        thetat = maxval;
        thetar = maxval1;

        ilt = maxnum;
        ilr = maxnum1;

        dlt = acc * ilt;
        dlr = d - acc * ilr;
    }
    else
    {
        thetat = thitai[pointnum];
        thetar = thitaj[0];
        ilt = (int)floor(((double)pointnum / 2.0));
        ilr = (int)ceil(((double)pointnum / 2.0));

        dlt = d / 2.0;  // dlt,dlr 未找到正确定义！
        dlr = d / 2.0;
    }
    theta = 1000 * d / ae + thetat + thetar;
    dct = d;  //！！！！待定
    dcr = d;

    //5.1.6 "光滑"地球模型和有效天线高度
    double Hst = 0;     // Hst: 路径起点，即干扰站上光滑地球表面的海拔高度
    double Hsr = 0;     // Hsr: 路径终点，即被干扰站上光滑地球表面的海拔高度

    Hst = ha - m * d / 2;
    Hsr = Hst + m * d / 2;

    if (Hst > hi[0])
    {
        Hst = hi[0];
        m = (Hsr - Hst) / d;
    }

    if (Hsr > hi[pointnum])
    {
        Hsr = hi[pointnum];
        m = (Hsr - Hst) / d;
    }

    Hte = Hts - Hst;
    Hre = Hrs - Hsr;

    //地形粗糙度hm
    double* hmi = new double[abs(ilr - ilt) + 1];
    hm = -1000;
    if (ilt <ilr)
    {
        for (int i = ilt; i < ilr + 1; i++)
        {
            hmi[i - ilt] = hi[i] - (Hst + m * di[i]);
            if (hm < hmi[i - ilt])
            {
                hm = hmi[i - ilt];
            }
        }
    }
    else
    {
        for (int i = ilr; i < ilt + 1; i++)
        {
            hmi[i - ilr] = hi[i] - (Hst + m * di[i]);
            if (hm < hmi[i - ilr])
            {
                hm = hmi[i - ilr];
            }
        }

    }
    Ag = AtmoGases(f, d, omiga);
	delete []Hci;
	delete []thitai;
	delete []thitaj;
	delete []hmi;
}

/*double attenuation_ra::final(double f, double theta, double p, double b0, double ae, double ab, double d, double N0, double Gt, double Gr, double omiga, double dlt, double dlr, double thetat, double thetar, int pointnum, double di[], double hi[])
{
    double Fj;  // 内插系数，考虑路径的角向距离
    double sita = 0.3;
    double sigma = 0.8;

    Fj = 1.0 - 0.5 * (1.0 + tanh(3.0 * sigma * (theta - sita) / sita));

    double Fk;  // 内插系数，考虑大圆路径长度
    double dsw = 20;    // 确定相关混合的距离范围的固定参数
    double K = 0.5;     // 确定在范围两端混合斜率的固定参数

    Fk = 1.0 - 0.5 * (1.0 + tanh(3.0 * K * (d - dsw) / dsw));

    // 计算与视距传播和海上部分路径绕射有关的假想最小传输损耗Lminb0p
    double L_minb0p;
    shiju(f, p, b0, d, dlt, dlr, Ag);

    raoshe(f, p, b0, omiga, d, pointnum, di, hi, ae, ab, Ag);

    if (p < b0)
        L_minb0p = L_b0p + (1 - omiga) * L_dp;
    else
        L_minb0p = L_bd50 + (L_b0b + (1 - omiga) * L_dp - L_bd50) * Fi;  // Fi 

    // 计算与视距和超视距信号增强有关的理论最小基本传输损耗 Lminbap
    double L_minbap;
    double ang = 2.5;
    double L_ba = atmoduct(f, omiga, p, b0, ae, d, dlt, dlr, theta, thetat, thetar, Ag);

    L_minbap = ang * log(exp(L_ba / ang) + exp(L_b0p / ang));

    // 计算与视距和超视距反射增强有关的理论基本传输损耗 Lbda
    double L_bda = (L_minbap > L_bd) ? L_bd : L_minbap + (L_bd - L_minbap) * Fk;

    // 计算修正的基本传输损耗 Lbam，该值纳入了绕射和视距或大气波导/层反射增强的影响
    double L_bam = L_bda + (L_minb0p - L_bda) * Fj;

    // 在 p% 时间内不超过的最终基本传输损耗Lb
    double L_bs = sanshe(f, p, d, theta, N0, Gt, Gr);
    double Aht = 0; // 散射体高度增益中引起的相应附加损耗
    double Ahr = 0;
    double L_b = -5 * log10(pow(10, -0.2 * L_bs) + pow(10, -0.2 * L_bam)) + Aht + Ahr;

    return L_b;
}*/

/*double attenuation_ra::Attenuation(double f, double omige, double p, double b0, double ae, double ab, double d, double di[], double hi[], int pointnum, double dii, double ha, double m, double N0, double Gt, double Gr, double x1, double y1, double x2, double y2,double h_radio)
{
    double L_final = 0; //总返回值
    double anzb1[2] = { x1, y1 };
    double anzb2[2] = { x2, y2 };
    if (x1 == x2 && y1 == y2)
        L_final = 0;

    else
    {
        maxnum = 0;
        maxnum1 = 0;
        maxval = -1000;
        maxval1 = -1000;
        dii = dii / 1000;

        paramcalcu(f, ae,d, di, hi, h_radio, pointnum, dii, ha, m, omige);

        L_final = final(f, theta, p, b0, ae, ab, d, N0, Gt, Gr, omige, dlt, dlr, thetat, thetar, pointnum, di, hi);


    }
    return L_final;
}*/
double chazhi(double** Dem, double czzb0, double czzb1, int Y, int X)
{
    double za = 0, zb = 0, zc = 0, zd = 0;
    double a = 0, b = 0, cz = 0;

    if (czzb0 >= Y || czzb1 >= X || czzb0 <= 1 || czzb1 <= 1) //当采样点为地图边界时，采用如下方法获得海拔高度值
    {
        if (czzb0 <= 1)
        {
            if (czzb0 <= 1 && czzb1 <= 1)
                cz = Dem[(int)(ceil(czzb0))][(int)(ceil(czzb1))];
            else if (czzb0 <= 1 && czzb1 >= X)
                cz = Dem[(int)(ceil(czzb0))][(int)(floor(czzb1)) - 2];
            else
                cz = Dem[(int)(ceil(czzb0))][(int)(ceil(czzb1)) - 1];
        }
        if (czzb1 <= 1)
        {
            if (czzb0 <= 1 && czzb1 <= 1)
                cz = Dem[(int)(ceil(czzb0))][(int)(ceil(czzb1))];
            else if (czzb0 >= Y && czzb1 <= 1)
                cz = Dem[(int)(floor(czzb0)) - 2][(int)(ceil(czzb1))];
            else
                cz = Dem[(int)(ceil(czzb0) - 1)][(int)(ceil(czzb1))];
        }
        if (czzb0 >= Y)
        {
            if (czzb0 >= Y && czzb1 <= 1)
                cz = Dem[(int)(floor(czzb0)) - 2][(int)(ceil(czzb1))];
            else if (czzb0 >= Y && czzb1 >= X)
                cz = Dem[(int)(floor(czzb0)) - 2][(int)(floor(czzb1)) - 2];
            else
                cz = Dem[(int)(floor(czzb0)) - 2][(int)(ceil(czzb1)) - 1];
        }
        if (czzb1 >= X)
        {
            if (czzb0 <= 1 && czzb1 >= X)
                cz = Dem[(int)(ceil(czzb0))][(int)(floor(czzb1)) - 2];
            else if (czzb0 >= Y && czzb1 >= X)
                cz = Dem[(int)(floor(czzb0)) - 2][(int)(floor(czzb1)) - 2];
            else
                cz = Dem[(int)(ceil(czzb0)) - 1][(int)(floor(czzb1)) - 2];
        }
    }

    else //当采样点为地图中间时，由以下内插方法计算得出海拔高度值
    {
        za = Dem[(int)floor(czzb0) - 1][(int)ceil(czzb1) - 1];
        zb = Dem[(int)ceil(czzb0) - 1][(int)ceil(czzb1) - 1];
        zc = Dem[(int)ceil(czzb0) - 1][(int)floor(czzb1) - 1];
        zd = Dem[(int)floor(czzb0) - 1][(int)floor(czzb1) - 1];
        a = czzb0 - floor(czzb0);
        b = czzb1 - floor(czzb1);
        cz = (1 - a) * b * za + a * b * zb + a * (1 - b) * zc + (1 - a) * (1 - b) * zd;
    }
    return cz;
}
double qixiang(double d, double fai)
{
    // fai：路径中心的纬度
    // 计算路径中心位置的异常传播的点发生率bita0
    double u1 = 0, u4 = 0.0;
    double tao = 0;
    double dtm = 0; // dtm为大圆路径的最长的连续的陆地（陆地加海滨）段（km）
    double dlm = 0;	// dlm为大圆路径的最长的连续的内陆段（km）P7表格
    double bita0;

    dtm = d;
    dlm = d;    // 信息不是很确定。
    tao = 1.0 - exp(-4.12 * pow(10.0, -4.0) * pow(dlm, 2.41));
    u1 = pow(pow(10.0, -dtm / (16.0 - 6.6 * tao)) + pow(pow(10.0, -(0.496 + 0.354 * tao)), 5.0), 0.2);
    if (fai <= 70)
    {	
        u4 = pow(10.0, (-0.935 + 0.0176 * fabs(fai)) * log10(u1));
        bita0 = pow(10.0, -0.015 * fabs(fai) + 1.67) * u1 * u4;
    }
    else
    {
        u4 = pow(10.0, 0.3 * log10(u1));
        bita0 = 4.17 * u1 * u4;
    }

    return bita0;
}
double pwtop(double pw, double fai, double w)
{
    // pw：最差月份的时间百分比
    // 时间百分数的计算
    double p;
    double Gl;

    if (fai <= 45)
        Gl = sqrt(1.1 + pow(fabs(cos(2 * fai * PI / 180)), 0.7));
    else
        Gl = sqrt(1.1 - pow(fabs(cos(2 * fai * PI / 180)), 0.7));

    p = pow(10, (log10(pw) + log10(Gl) - 0.186 * w - 0.444) / (0.816 + 0.078 * w));
    return p;
}
double attenuation_ra::Radio_C(double** Dem, double deltaN, double xt, double yt, double xr, double yr, int rows, int columns,double fre,double Gt,double Gr, double acc,double dx,double dy, double pw,double fai,double height1)
{
	double L_final; //总返回值
    // pw：最差月份的时间百分比 //  fai：路径中心的纬度
	double omige =0.0;
	double ae; // 地球等效半径
	double ab; // 地球等效半径
	double d;  // 通信路径的距离
	int pointnum; // 通信链路的取样点数
	double *di; // 从发射端到每个取样点的距离
	double *hi; // 每个取样点的高程值

	double b0 = 0; // 异常传播发生的概率
	double p = 0; // 选择的时间百分比（平均年份预测中使用的时间百分比）
	double ha = 0; // 实际传播路径海拔高度的平均值
	double m = 0; // 相对于海平面的最小二乘方表面的斜率

	lamde=0;		// 微波的波长（km）
    h1=0;			// 发射天线的高度
    h2=0;			// 计算出来的接收天线的高度
    Hts=0;			// 发射端天线的平均高度
    Hrs=0;			// 接收端天线的平均高度
    maxval=0;		// 发射端的物理视界仰角
    maxnum=0;		// 发射端视界临界点在路径上的位置
    maxval1=0;		// 接收端的物理视界仰角
    maxnum1=0;		// 接收端视界临界点在路径上的位置
    type=0;			// 视界or非视界的标志位

    thetat=0, thetar=0;	// 发射端（接收端）视界俯仰角
    theta=0;			// 路径的角向距离
    ilt=0, ilr=0;			// 发射端（接收端）视界临界点在链路上的位置
    dlt=0 , dlr=0;		// 从发射天线和接收天线到相应的视界的距离。
	dct=0 , dcr=0;		// 从发射端（接收端）到沿大圆路径上跨越陆地的距离

    Hte=0;   // 发射端天线有效高度
    Hre=0;   // 接收端天线的有效高度
    hm =0;   // 地面粗糙度
    L_b0p=0;		// 视距传播P%时间内不超过的基本损耗
	L_b0b=0;		// 视距传播b0%时间内不超过的基本损耗
    Ag=0;			// 大气吸收损耗

        //raoshe raoshe raoshe raoshe raoshe raoshe raoshe raoshe raoshe raoshe raoshe raoshe raoshe raoshe raoshe 
    Ld50=0;  // 绕射损耗中值Ld50√
    Ldb=0;       //b%时间内不超过的绕射损耗Ldβ √
    L_dp=0;      //对于p%时间不超过的绕射损耗√
    L_bd50=0;    //与绕射有关的基本传输损耗中值√
    L_bd=0;      //与绕射有关的对于p%时间不超过的基本传输损耗√
    Fi=0;    //√

    if ((xt == xr) && (yt == yr))
    {
        ae = 0;
        ab = 0;
        m = 0;
        p = 0;
        pointnum = 0;
        b0 = 0;
        hi = new double[0];
        ha = 0;
        di = new double[0];
        d = 0;
    }
    else
    {
        //double dx = 73.5624;//double dy = 92.6043;	
        double a = 6371;		//地球半径 km
        double K50 = 157 / (157 - deltaN);  // 有效地球半径因子中值
        double Kb = 3.0;    // bita0%时间内超过的有效地球半径系数估计值

		//double fai=37.0;//( yt + yr ) / 2.0;
		int i;
        //中间参数，计算用
        double dzbx = 0, dzby = 0;
        //===========有效地球半径计算===================
        ae = a * K50;
        ab = a * Kb;

        d = sqrt((xt - xr) * dx*(xt - xr) * dx + (yt - yr) * dy*(yt - yr) * dy) / 1000.0;

        b0 = qixiang(d, fai);
        p = pwtop(pw, fai, 0);

        pointnum = (int)(d / acc);
		if (pointnum == 0)
        {
            pointnum = 1;
        }

        dzbx = (xt - xr) / pointnum;
        dzby = (yt - yr) / pointnum;

        di = new double[pointnum + 1];
		double **czzb;
		czzb = new double*[pointnum + 1];
		for(i=0;i < pointnum + 1;i++)
		{
			czzb[i] = new double[2]; 
		}
        hi = new double[pointnum + 1];

        for (i = 0; i < pointnum + 1; i++)  //√
        {
            di[i] = i * acc;
            czzb[i][0] = xt - i * dzbx;
            czzb[i][1] = yt - i * dzby;
            hi[i] = chazhi(Dem, czzb[i][0], czzb[i][1],rows, columns);
        }
		for(i=0; i<pointnum + 1;i++)
		{
			delete[]czzb[i];
		}
		delete []czzb;
        double sum1 = 0, sum2 = 0,sum3 = 0;  // 计算m用
        for (i = 0; i < pointnum + 1; i++)
        {
            sum3 += hi[i];
        }
        ha = sum3 / (pointnum + 1);

        for (i = 0; i < pointnum + 1; i++)
        {
            sum1 += (hi[i] - ha) * (di[i] - d / 2);
            sum2 += (di[i] - d / 2) * (di[i] - d / 2);
        }
        m = sum1 / sum2;
    }

	if ((xt == xr) && (yt == yr))
	{
		L_final=0.0;
	}
	else
	{
        maxnum = 0;
        maxnum1 = 0;
        maxval = -1000;
        maxval1 = -1000;
		paramcalcu(fre, ae,d, di, hi, height1, pointnum, acc, ha, m, omige);
	// ========================== 内插系数的计算 ===================================== //
		double Fj;          // 内插系数，考虑路径的角向距离
		double sita = 0.3;
		double sigma = 0.8;
		double Fk;          // 内插系数，考虑大圆路径长度
		double dsw = 20;    // 确定相关混合的距离范围的固定参数
		double K = 0.5;     // 确定在范围两端混合斜率的固定参数
		Fj = 1.0 - 0.5 * (1.0 + tanh(3.0 * sigma * (theta - sita) / sita));
		Fk = 1.0 - 0.5 * (1.0 + tanh(3.0 * K * (d - dsw) / dsw));
		double Aht = 0.0; // 散射体高度增益中引起的相应附加损耗
		double Ahr = 0.0; // 散射体高度增益中引起的相应附加损耗
	// ========================== 传输损耗预测算 ===================================== //
		if(type == 0)
		{
			shiju(fre, p, b0, d, dlt, dlr, Ag);
			L_final = L_b0p + Aht + Ahr;
		}
		else
		{
			shiju(fre, p, b0, d, dlt, dlr, Ag);
			raoshe(fre, p, b0, omige, d, pointnum, di, hi, ae, ab, Ag);

			// 计算与视距传播和海上部分路径绕射有关的假想最小传输损耗Lminb0p
			double L_minb0p;
			if (p < b0)
				L_minb0p = L_b0p + (1 - omige) * L_dp;
			else
				L_minb0p = L_bd50 + (L_b0b + (1 - omige) * L_dp - L_bd50) * Fi;  // Fi 

			// 计算与视距和超视距信号增强有关的理论最小基本传输损耗 Lminbap
			double L_minbap;
			double ang = 2.5;
			double L_ba = atmoduct(fre, omige, p, b0, ae, d, dlt, dlr, theta, thetat, thetar, Ag);
			L_minbap = ang * log(exp(L_ba / ang) + exp(L_b0p / ang));

			// 计算与视距和超视距反射增强有关的理论基本传输损耗 Lbda
			double L_bda = (L_minbap > L_bd) ? L_bd : L_minbap + (L_bd - L_minbap) * Fk;

			// 计算修正的基本传输损耗 Lbam，该值纳入了绕射和视距或大气波导/层反射增强的影响
			double L_bam = L_bda + (L_minb0p - L_bda) * Fj;

			// 在 p% 时间内不超过的最终基本传输损耗L_final
			double L_bs = sanshe(fre, p, d, theta, 325.0, Gt, Gr);
			L_final = -5 * log10(pow(10, -0.2 * L_bs) + pow(10, -0.2 * L_bam)) + Aht + Ahr;

		}	

	}	
	delete []di;
	delete []hi;

	return L_final;
}