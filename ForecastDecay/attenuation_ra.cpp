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
static double lamde;		// ΢���Ĳ�����km��
static double h1;			// �������ߵĸ߶�
static double h2;			// ��������Ľ������ߵĸ߶�
static double Hts;			// ��������ߵ�ƽ���߶�
static double Hrs;			// ���ն����ߵ�ƽ���߶�
static double maxval;		// ����˵������ӽ�����
static int maxnum;			// ������ӽ��ٽ����·���ϵ�λ��
static double maxval1;		// ���ն˵������ӽ�����
static int maxnum1;			// ���ն��ӽ��ٽ����·���ϵ�λ��
static int type;			// �ӽ�or���ӽ�ı�־λ

static double thetat, thetar;	// ����ˣ����նˣ��ӽ縩����
static double theta;			// ·���Ľ������
static int ilt, ilr;			// ����ˣ����նˣ��ӽ��ٽ������·�ϵ�λ��
static double dlt , dlr;		// �ӷ������ߺͽ������ߵ���Ӧ���ӽ�ľ��롣
static double dct , dcr;		// �ӷ���ˣ����նˣ����ش�Բ·���Ͽ�Խ½�صľ���

static double Hte;				// �����������Ч�߶�
static double Hre;				// ���ն����ߵ���Ч�߶�
static double hm ;				// ����ֲڶ�
// = = = = = = = = = = = = = = = �Ӿ� = = = = = = = = = = = = = = = = = = = = = = = = //  
static double L_b0p;			// �Ӿഫ��P%ʱ���ڲ������Ļ������
static double L_b0b;			// �Ӿഫ��b0%ʱ���ڲ������Ļ������
static double Ag;				// �����������

// = = = = = = = = = = = = = = = ���� = = = = = = = = = = = = = = = = = = = = = = = = //
static double Ld50;				// ���������ֵLd50��
static double Ldb;				// b%ʱ���ڲ��������������Ld�� ��
static double L_dp;				// ����p%ʱ�䲻������������ġ�
static  double L_bd50;			// �������йصĻ������������ֵ��
static  double L_bd;			// �������йصĶ���p%ʱ�䲻�����Ļ���������ġ�
static double Fi;				// ��
void attenuation_ra::shiju(double f,double p,double b0, double d, double dlt, double dlr, double Ag)
{
	double L_bfsg; // �Ӿഫ���Ļ����������
	double Esp,Esb;
	L_bfsg = 92.5 + 20 * log10(f) + 20 *log10(d) +Ag;
	Esp = 2.6 * (1.0 - exp(-0.1*(dlt+dlr))) * log10(p/50.0); // ��p%ʱ�����ɶྶ�;۽�ЧӦ�����У����
	Esb = 2.6 * (1.0 - exp(-0.1*(dlt+dlr))) * log10(b0/50.0); // ��b0%ʱ�����ɶྶ�;۽�ЧӦ�����У����
	L_b0p = L_bfsg + Esp; // ��p%ʱ���ڲ����������Ӿഫ������Ļ����������
	L_b0b = L_bfsg + Esb; // ��b0%ʱ���ڲ����������Ӿഫ������Ļ����������
}

double calcuI(double x)
{
	// ���ۻ���̫����������������ĳɶ�����̫�ֲ����ڲ�����
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
	// ��Ш��������Ľ���ֵ�ĺ���
    double J;
    if (v > -0.78)
        J = 6.9 + 20 * log10(sqrt((v - 0.1) * (v - 0.1) + 1) + v - 0.1);
    else J = 0;
    return J;
}

void attenuation_ra::raosheParam(int pointnum, double di[], double hi[], double d, double ae, double abeite)
{
    double keseim = 0;  //��������·��б�ʵ�У�����m
    double* Hi = new double[pointnum + 1];  //��ֱ����Hi
    double* vmi50 = new double[pointnum + 1];
    double vm50 = -1000;    //�������
    int im50 = 0;   //�������ֵvm50��������ָ��

    keseim = cos(atan(1e-3 * (Hrs - Hts) / d));
    // ============================ ��������·�����������߼����С���============================ //
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
	// ================================ �������������ֵ ========================================== //
    double Lm50 = 0;    // ��Ҫ�߽��Ш�����������ֵLm50
    double Lt50 = 0;    // �������α߽��Ш�����������ֵLt50
    double Lr50 = 0;    // ���ջ���α߽��Ш�����������ֵLr50

    double keseit = cos(atan(1e-3 * (Hi[im50] - Hts) / di[im50]));          // У����
    double keseir = cos(atan(1e-3 * (Hrs - Hi[im50]) / (d - di[im50])));    // У����
    int it50 = 0;   // �������α߽�������ָ��
    int ir50 = 0;   // ���ջ���α߽�������ָ��

    Lm50 = calcuJ(vm50);
    // ��Lm50 = 0�������������ֵLd50��Ld�¾�Ϊ0�������ټ������䡣
    if (Lm50 == 0)
    {
        Ld50 = 0;
        Ldb = 0;
    }
    // ��Lm50 �� 0�������¼������߽�ķ�����ͽ��ջ���Ĵα߽����������������
    else
    {
        //------------------------- �������α߽�����������ֵLt50 ----------------------------- //
        if (im50 == 1)   // �����ڷ������α߽磬����������Lt50����
            Lt50 = 0;
        else
        {
            double* vti50 = new double[pointnum + 1];
            double vt50 = -1000;    // �������

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
        //----------------------------- ���ջ���α߽�����������ֵLr50 ----------------------------- //
        if (im50 == pointnum - 1)
            Lr50 = 0;
        else
        {
            double* vri50 = new double[pointnum + 1];
            double vr50 = -1000;    //�������

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

        // 4.2.2 �ڦ�0% ʱ���ڲ��������������Ld��
        // ������Ч����뾶 a��������
        // �ڦ�0%ʱ���ڲ����������߽�Ш���������Lm��
        double Lmbita;
        double vmbita;  // �������vm��
        double Himbita = hi[im50] + 1e3 * di[im50] * (d - di[im50]) / (2 * abeite) - (Hts * (d - di[im50]) + Hrs * di[im50]) / d;

        vmbita = keseim * Himbita * sqrt(2e-3 * d / (lamde * di[im50] * (d - di[im50])));
        Lmbita = calcuJ(vmbita);

        // �ڦ�0%ʱ���ڲ������ķ������α߽�Ш���������Lt��
        double Ltbita;

        // �� Lt50 = 0��Lt�½�Ϊ0����֮������Lt�����£�
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

        // �ڦ�0%ʱ���ڲ������Ľ��ջ���α߽�Ш���������Lr��
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
        // �ڦ�0%ʱ���ڲ������ĸ��ֱ߽���������һ��
        Ldb = Lmbita + (1 - exp(-Lmbita / 6)) * (Ltbita + Lrbita + 10 + 0.04 * d);
    }
	
	
	delete []vmi50;
	delete []Hi;
}

void attenuation_ra::raoshe(double f, double p, double beite0, double omiga, double d, int pointnum, double di[], double hi[], double ae, double abeite, double Ag)
{
    // ������ģ�ͼ����������Ҫ����3��ֵ��
    // ����p%ʱ�䲻������������� L_dp
    // �������йصĻ������������ֵ L_bd50
    // �������йصĶ���p%ʱ�䲻�����Ļ���������� L_bd 

    double L_bfsg = 0;
    raosheParam(pointnum, di, hi, d, ae, abeite);

    L_bfsg = 92.5 + 20 * log10(f) + 20 * log10(d) + Ag;

    if (p == 50)
        Fi = 0;
    else if (p < 50 && p > beite0)
        Fi = calcuI(p / 100) / calcuI(beite0 / 100);
    else
        Fi = 1;

    L_dp = Ld50 + Fi * (Ldb - Ld50); // ����p%ʱ�䲻������������� L_dp
    L_bd50 = L_bfsg + Ld50;          // �������йصĻ������������ֵ L_bd50
    L_bd = L_b0p + L_dp;             // �������йصĶ���p%ʱ�䲻�����Ļ���������� L_bd 
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
    // �Ӻ�ƽ�浽10 km�߶ȷ�Χ�ڣ��ɸ��������ˮ����ɵ����ߵ粨����˥��
    double Ag = 0;  // �ܵ��������գ�����ֵ
    double ro = 0;  // ����˥����o(dB/km)
    double rw = 0;  // ����ˮ���е�˥��������˥��ֵ��w (dB/km)

    // �м����
    double rt = 0, rp = 0;
    double kesei1 = 0, kesei2 = 0, kesei3 = 0;  // ���� ro �м����
    double yita1 = 0, yita2 = 0;  // ���� rw �м����

    double pres = 1013.25; // ��ѹ��hPa��
    double temp = 15;   // �¶ȣ���C��
    double rou = 7.5 + 2.5 * omiga;    // ˮ�����ܶ�(g/m^3)

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
    double Lf, Lc, Ag = 0;	// AgΪ��������

    Ag = AtmoGases(f, d, -1.8);
    Lf = 25 * log10(f) - 2.5 * log10(f / 2.0) * log10(f / 2.0);	// ��Ƶ���йص����
    Lc = 0.051 * exp(0.055 * (Gt + Gr));
    L_bs = 190 + Lf + 20 * log10(d) + 0.573 * thita - 0.15 * N0 + Lc + Ag - 10.1 * pow((-log10(p / 50)), 0.7);
    return L_bs;
}

double attenuation_ra::atmoduct(double f, double omiga, double p, double b0, double ae, double d, double dlt, double dlr, double thita, double thitat, double thitar, double Ag)
{
    double L_ba;	// ��������/�㷴�����
    double A_f;     // �����������ߺ��쳣�����ṹ֮��Ĺ̶������ĵ��ۺϣ�����ɢ����ĳ��⣩
    double A_dp;	// ���쳣������������ʱ��ٷֱȺͽǶ�--�����йص����

    // ===============================================1. ����A_f================================================== //
    // �м����
    double Ast = 0, Asr = 0;	// ����վ��������վ��λ���������
    double Act = 0, Acr = 0;	// ����վ��������վ�Ŀ纣��������������У����
    double thitatpp, thitarpp;

    thitatpp = thitat - 0.1 * dlt;    // mrad
    thitarpp = thitar - 0.1 * dlr;    // mrad

    // ------------------------------------------ ����Ast��Asr ----------------------------------------- //
    if (thitatpp <= 0)
        Ast = 0;
    else
        Ast = 20 * log10(1 + 0.361 * thitatpp * sqrt(f * dlt)) + 0.264 * thitatpp * pow(f, 1 / 3);

    if (thitarpp <= 0)
        Asr = 0;
    else
        Asr = 20 * log10(1 + 0.361 * thitarpp * sqrt(f * dlr)) + 0.264 * thitarpp * pow(f, 1 / 3);

    // ------------------------------------------ ����Act��Acr ------------------------------------------ //
    if (omiga >= 0.75 && dct <= dlt && dct <= 5)
        Act = -3 * exp(-0.25 * dct * dct) * (1 + tanh(0.07 * (50 - Hts)));
    else
        Act = 0;

    if (omiga >= 0.75 && dcr <= dlr && dcr <= 5)
        Acr = -3 * exp(-0.25 * dcr * dcr) * (1 + tanh(0.07 * (50 - Hrs)));
    else
        Acr = 0;

    A_f = 102.45 + 20 * log10(f) + 20 * log10(dlt + dlr) + Ast + Asr + Act + Acr; //!!!! dlt,dlr 

    //============================================2. ����A_dp===============================================//
    // ------------------------------------------ 2.1 ����rd ------------------------------------------ //
    double rd;	// ��˥��
    rd = 5e-5 * ae * pow(f, 1 / 3);

    // ------------------------------------------ 2.2 ���� thitap -------------------------------------- //
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

    // ------------------------------------------ 2.3 ����Ap ------------------------------------------ //
    double Ap;	//ʱ��ٷֱ�p�Ŀɱ��ԣ��ۻ��ֲ���
    double TT;  // TT : ��
    double bb;
    double u2, u3;
    double aa;
    double dd = 3.5;   // dd:��
    double dlm = 0;    // ��Բ·�����������½�ء�����ȫ��½�ص������dlm = d
    double tal = floor(1 - exp(-4.12e-4 * pow(dlm, 2.41)));    //��ʽ 3-133a���Ƿ���ȡ��������
    double dI = (d - dlt - dlr < 40) ? d - dlt - dlr : 40;
    double Temp = 0.0;

    aa = -0.6 - dd * 1e-9 * pow(d, 3.1) * tal;
    Temp = 500 * d * d / (ae * pow(sqrt( Hte) + sqrt( Hre), 2));
    u2 = pow(500 * d * d / (ae * pow(sqrt( Hte) + sqrt( Hre), 2)), aa);   //!!! ��ϵ����Ҫ����

    u2 = (u2 > 1) ? 1 : u2;

    if ( hm <= 10)
        u3 = 1;
    else
        u3 = exp(-4.6e-5 * ( hm - 10) * (43 + 6 * dI));

    bb = b0 * u2 * u3;

    TT = 1.076 / pow(2.0058 - log10(bb), 1.012) * exp(-(9.51 - 4.8 * log10(bb) + 0.198 * log10(bb) * log10(bb)) * 1e-6 * pow(d, 1.13));
    Ap = -12 + (1.2 + 3.7e-3 * d) * log10(p / bb) + 12 * pow(p / bb, TT);
	// �쳣������������ʱ��ٷ����ͽǶȡ������йص����
    A_dp = rd * thitap + Ap;
	// �쳣�����ڼ�Ļ����������
    L_ba = A_f + A_dp + Ag;

    return L_ba;
}
  
void attenuation_ra::paramcalcu(double f, double ae, double d, double di[], double hi[], double h_radio, int pointnum, double acc,double ha, double m, double omiga)
{
    lamde = 0.3 / f;
	int i = 0;
    int n = -1;
    int hight = (int)floor(h_radio);
    h1 = 10;	// �������߸߶�
    h2 = 0;		// �������߸߶�
    double* Hci = new double[pointnum + 1];    //��ֵ����϶
    double Hcmin = -100;   //��С��϶
   // ==================================== ���߸߶ȼ��� ========================================== //
	// ���㴫����϶
    for (i = 1; i < pointnum; i++)
    {
        Hci[i] = ((hi[0] + hight) * (d - di[i]) + (hi[pointnum] + hight) * di[i]) / d - hi[i];
    }
	 double minHci = Hci[1];	// ���µ�·����϶
    int minHcNum = 1;			// ��С��϶���λ��
    double minHr = 0;			// �������ߵ���С�߶�
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

	// ==================================== ���߸߶ȼ������ ========================================== //
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
        thitai[i] = (hi[i] - Hts) / di[i] - 1000 * di[i] / (2 * ae);	//thitai ��λΪ mrad
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

        dlt = d / 2.0;  // dlt,dlr δ�ҵ���ȷ���壡
        dlr = d / 2.0;
    }
    theta = 1000 * d / ae + thetat + thetar;
    dct = d;  //������������
    dcr = d;

    //5.1.6 "�⻬"����ģ�ͺ���Ч���߸߶�
    double Hst = 0;     // Hst: ·����㣬������վ�Ϲ⻬�������ĺ��θ߶�
    double Hsr = 0;     // Hsr: ·���յ㣬��������վ�Ϲ⻬�������ĺ��θ߶�

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

    //���δֲڶ�hm
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
    double Fj;  // �ڲ�ϵ��������·���Ľ������
    double sita = 0.3;
    double sigma = 0.8;

    Fj = 1.0 - 0.5 * (1.0 + tanh(3.0 * sigma * (theta - sita) / sita));

    double Fk;  // �ڲ�ϵ�������Ǵ�Բ·������
    double dsw = 20;    // ȷ����ػ�ϵľ��뷶Χ�Ĺ̶�����
    double K = 0.5;     // ȷ���ڷ�Χ���˻��б�ʵĹ̶�����

    Fk = 1.0 - 0.5 * (1.0 + tanh(3.0 * K * (d - dsw) / dsw));

    // �������Ӿഫ���ͺ��ϲ���·�������йصļ�����С�������Lminb0p
    double L_minb0p;
    shiju(f, p, b0, d, dlt, dlr, Ag);

    raoshe(f, p, b0, omiga, d, pointnum, di, hi, ae, ab, Ag);

    if (p < b0)
        L_minb0p = L_b0p + (1 - omiga) * L_dp;
    else
        L_minb0p = L_bd50 + (L_b0b + (1 - omiga) * L_dp - L_bd50) * Fi;  // Fi 

    // �������Ӿ�ͳ��Ӿ��ź���ǿ�йص�������С����������� Lminbap
    double L_minbap;
    double ang = 2.5;
    double L_ba = atmoduct(f, omiga, p, b0, ae, d, dlt, dlr, theta, thetat, thetar, Ag);

    L_minbap = ang * log(exp(L_ba / ang) + exp(L_b0p / ang));

    // �������Ӿ�ͳ��Ӿ෴����ǿ�йص����ۻ���������� Lbda
    double L_bda = (L_minbap > L_bd) ? L_bd : L_minbap + (L_bd - L_minbap) * Fk;

    // ���������Ļ���������� Lbam����ֵ������������Ӿ���������/�㷴����ǿ��Ӱ��
    double L_bam = L_bda + (L_minb0p - L_bda) * Fj;

    // �� p% ʱ���ڲ����������ջ����������Lb
    double L_bs = sanshe(f, p, d, theta, N0, Gt, Gr);
    double Aht = 0; // ɢ����߶��������������Ӧ�������
    double Ahr = 0;
    double L_b = -5 * log10(pow(10, -0.2 * L_bs) + pow(10, -0.2 * L_bam)) + Aht + Ahr;

    return L_b;
}*/

/*double attenuation_ra::Attenuation(double f, double omige, double p, double b0, double ae, double ab, double d, double di[], double hi[], int pointnum, double dii, double ha, double m, double N0, double Gt, double Gr, double x1, double y1, double x2, double y2,double h_radio)
{
    double L_final = 0; //�ܷ���ֵ
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

    if (czzb0 >= Y || czzb1 >= X || czzb0 <= 1 || czzb1 <= 1) //��������Ϊ��ͼ�߽�ʱ���������·�����ú��θ߶�ֵ
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

    else //��������Ϊ��ͼ�м�ʱ���������ڲ巽������ó����θ߶�ֵ
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
    // fai��·�����ĵ�γ��
    // ����·������λ�õ��쳣�����ĵ㷢����bita0
    double u1 = 0, u4 = 0.0;
    double tao = 0;
    double dtm = 0; // dtmΪ��Բ·�������������½�أ�½�ؼӺ������Σ�km��
    double dlm = 0;	// dlmΪ��Բ·���������������½�Σ�km��P7���
    double bita0;

    dtm = d;
    dlm = d;    // ��Ϣ���Ǻ�ȷ����
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
    // pw������·ݵ�ʱ��ٷֱ�
    // ʱ��ٷ����ļ���
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
	double L_final; //�ܷ���ֵ
    // pw������·ݵ�ʱ��ٷֱ� //  fai��·�����ĵ�γ��
	double omige =0.0;
	double ae; // �����Ч�뾶
	double ab; // �����Ч�뾶
	double d;  // ͨ��·���ľ���
	int pointnum; // ͨ����·��ȡ������
	double *di; // �ӷ���˵�ÿ��ȡ����ľ���
	double *hi; // ÿ��ȡ����ĸ߳�ֵ

	double b0 = 0; // �쳣���������ĸ���
	double p = 0; // ѡ���ʱ��ٷֱȣ�ƽ�����Ԥ����ʹ�õ�ʱ��ٷֱȣ�
	double ha = 0; // ʵ�ʴ���·�����θ߶ȵ�ƽ��ֵ
	double m = 0; // ����ں�ƽ�����С���˷������б��

	lamde=0;		// ΢���Ĳ�����km��
    h1=0;			// �������ߵĸ߶�
    h2=0;			// ��������Ľ������ߵĸ߶�
    Hts=0;			// ��������ߵ�ƽ���߶�
    Hrs=0;			// ���ն����ߵ�ƽ���߶�
    maxval=0;		// ����˵������ӽ�����
    maxnum=0;		// ������ӽ��ٽ����·���ϵ�λ��
    maxval1=0;		// ���ն˵������ӽ�����
    maxnum1=0;		// ���ն��ӽ��ٽ����·���ϵ�λ��
    type=0;			// �ӽ�or���ӽ�ı�־λ

    thetat=0, thetar=0;	// ����ˣ����նˣ��ӽ縩����
    theta=0;			// ·���Ľ������
    ilt=0, ilr=0;			// ����ˣ����նˣ��ӽ��ٽ������·�ϵ�λ��
    dlt=0 , dlr=0;		// �ӷ������ߺͽ������ߵ���Ӧ���ӽ�ľ��롣
	dct=0 , dcr=0;		// �ӷ���ˣ����նˣ����ش�Բ·���Ͽ�Խ½�صľ���

    Hte=0;   // �����������Ч�߶�
    Hre=0;   // ���ն����ߵ���Ч�߶�
    hm =0;   // ����ֲڶ�
    L_b0p=0;		// �Ӿഫ��P%ʱ���ڲ������Ļ������
	L_b0b=0;		// �Ӿഫ��b0%ʱ���ڲ������Ļ������
    Ag=0;			// �����������

        //raoshe raoshe raoshe raoshe raoshe raoshe raoshe raoshe raoshe raoshe raoshe raoshe raoshe raoshe raoshe 
    Ld50=0;  // ���������ֵLd50��
    Ldb=0;       //b%ʱ���ڲ��������������Ld�� ��
    L_dp=0;      //����p%ʱ�䲻������������ġ�
    L_bd50=0;    //�������йصĻ������������ֵ��
    L_bd=0;      //�������йصĶ���p%ʱ�䲻�����Ļ���������ġ�
    Fi=0;    //��

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
        double a = 6371;		//����뾶 km
        double K50 = 157 / (157 - deltaN);  // ��Ч����뾶������ֵ
        double Kb = 3.0;    // bita0%ʱ���ڳ�������Ч����뾶ϵ������ֵ

		//double fai=37.0;//( yt + yr ) / 2.0;
		int i;
        //�м������������
        double dzbx = 0, dzby = 0;
        //===========��Ч����뾶����===================
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

        for (i = 0; i < pointnum + 1; i++)  //��
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
        double sum1 = 0, sum2 = 0,sum3 = 0;  // ����m��
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
	// ========================== �ڲ�ϵ���ļ��� ===================================== //
		double Fj;          // �ڲ�ϵ��������·���Ľ������
		double sita = 0.3;
		double sigma = 0.8;
		double Fk;          // �ڲ�ϵ�������Ǵ�Բ·������
		double dsw = 20;    // ȷ����ػ�ϵľ��뷶Χ�Ĺ̶�����
		double K = 0.5;     // ȷ���ڷ�Χ���˻��б�ʵĹ̶�����
		Fj = 1.0 - 0.5 * (1.0 + tanh(3.0 * sigma * (theta - sita) / sita));
		Fk = 1.0 - 0.5 * (1.0 + tanh(3.0 * K * (d - dsw) / dsw));
		double Aht = 0.0; // ɢ����߶��������������Ӧ�������
		double Ahr = 0.0; // ɢ����߶��������������Ӧ�������
	// ========================== �������Ԥ���� ===================================== //
		if(type == 0)
		{
			shiju(fre, p, b0, d, dlt, dlr, Ag);
			L_final = L_b0p + Aht + Ahr;
		}
		else
		{
			shiju(fre, p, b0, d, dlt, dlr, Ag);
			raoshe(fre, p, b0, omige, d, pointnum, di, hi, ae, ab, Ag);

			// �������Ӿഫ���ͺ��ϲ���·�������йصļ�����С�������Lminb0p
			double L_minb0p;
			if (p < b0)
				L_minb0p = L_b0p + (1 - omige) * L_dp;
			else
				L_minb0p = L_bd50 + (L_b0b + (1 - omige) * L_dp - L_bd50) * Fi;  // Fi 

			// �������Ӿ�ͳ��Ӿ��ź���ǿ�йص�������С����������� Lminbap
			double L_minbap;
			double ang = 2.5;
			double L_ba = atmoduct(fre, omige, p, b0, ae, d, dlt, dlr, theta, thetat, thetar, Ag);
			L_minbap = ang * log(exp(L_ba / ang) + exp(L_b0p / ang));

			// �������Ӿ�ͳ��Ӿ෴����ǿ�йص����ۻ���������� Lbda
			double L_bda = (L_minbap > L_bd) ? L_bd : L_minbap + (L_bd - L_minbap) * Fk;

			// ���������Ļ���������� Lbam����ֵ������������Ӿ���������/�㷴����ǿ��Ӱ��
			double L_bam = L_bda + (L_minb0p - L_bda) * Fj;

			// �� p% ʱ���ڲ����������ջ����������L_final
			double L_bs = sanshe(fre, p, d, theta, 325.0, Gt, Gr);
			L_final = -5 * log10(pow(10, -0.2 * L_bs) + pow(10, -0.2 * L_bam)) + Aht + Ahr;

		}	

	}	
	delete []di;
	delete []hi;

	return L_final;
}