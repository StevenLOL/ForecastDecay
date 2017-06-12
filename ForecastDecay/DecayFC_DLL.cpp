#include "DecayFC_DLL.h"
#include <stdio.h>
#include <iostream>
#include <math.h>
#include "ITUrp.h"
#include "attenuation_ra.h"
using namespace std;
#define PI 3.14159265358

double **DecayDiffraction(double Setparameter[], double Receparameter[], Data **Eledata, int dimension[], double f, Iparameter IPara, double dx1, double dy1, int mark)
{

	double xt1, yt1, Gt, ht, xr1, yr1, Gr, hr, acc, aa, distance;
	int row, column, e, i, j;
	xt1 = Setparameter[0]; // �������߾���
	yt1 = Setparameter[1]; // ��������γ��
	Gt = Setparameter[2]; // �������������������
	ht = Setparameter[3]; // �������߸߶�

	xr1 = Receparameter[0]; // �������߾���
	yr1 = Receparameter[1]; // ��������γ��
	Gr = Receparameter[2]; // �������������������
	hr = Receparameter[3]; // �������߸߶�

	row = dimension[0];    // ������Ϣ���ݵ�����
	column = dimension[1]; // ������Ϣ���ݵ�����
	acc = IPara.pSize;	   // Ԥ�������
	distance = IPara.Elength;// ��ЧԤ�ⷶΧ
	if (IPara.aMark == 0)// ��ˮ
	{
		aa = 0.001;		// ��Ч�絼��
		e = 80;			// ��Խ�糣��
	}
	if (IPara.aMark == 1)// ��ˮ
	{
		aa = 4.0;			// ��Ч�絼��
		e = 80;			// ��Խ�糣��
	}
	if (IPara.aMark == 2)// ʪ��
	{
		aa = 0.008;		// ��Ч�絼��
		e = 10;			// ��Խ�糣��
	}
	if (IPara.aMark == 3)// �ɵ�
	{
		aa = 0.001;		// ��Ч�絼��
		e = 2;			// ��Խ�糣��
	}
	double **Fvalue; // ���巵��ֵ��άָ��
	double **GRID;   // ����߳̾���ָ��
	int y3 = 1, y4 = 0;		// ���������͵������
	double re = 8500000;  // �����Ч�뾶
	// ======================== ��ͼ�����ڲ������ˮƽ���ʹ�ֱ���m ============================ //

	// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = //
	GRID = new double*[row];
	for (i = 0; i < row; i++)
	{
		GRID[i] = new double[column];
	}
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < column; j++)
		{
			GRID[i][j] = Eledata[i][j].eleva;
		}
	}
	for (i = 0; i<row; i++)
	{
		delete[]Eledata[i];
	}
	delete[]Eledata;

 //mark=0  ����ǿԤ��
	if (mark == 0)
	{
		Fvalue = new double*[row];
			for (i = 0; i < row; i++)
			{
				Fvalue[i] = new double[column];
			}
			//////////////////////// ITU-rp526��ǿԤ��ģ���� ////////////////////////////////////////////


			ITUrp RaosheITU;
			double lined;
			for (i = 0; i < row; i++)
			{
				for (j = 0; j < column; j++)
				{
					lined = sqrt(((xt1 - i) * dx1) * ((xt1 - i) * dx1) + ((yt1 - j) * dy1) * ((yt1 - j) * dy1));//����վ����ͼĳ��ľ���
					if (lined == 0)
					{
						Fvalue[i][j] = 0.0;
					}
					//else if (lined > distance * 1000)// ��������Ϊ����ͨ��
					//{
					//	Fvalue[i][j] = 1000;	// ����ͨ�ŷ�Χ����Ϊ��ļ�����Ϊ˥��1000db
					//}
					else if (lined < 1000)
					{
						Fvalue[i][j] = 32.45 + 20 * log10(lined / 1000) + 26 * log10(f); //���ɿռ䴫��
					}
					else
					{
						Fvalue[i][j] = RaosheITU.ITUraoshe(xt1, yt1, i, j, GRID, row, column, acc, f, ht, hr, y3, y4, aa, e, re, dx1, dy1);
					}
				}

			}
	}
	else
	{
		//mark=1 ��·˥��Ԥ��
		double Length_L, lined, spacX, spacY, xi, yj;
		int N_Sample;
		Length_L = sqrt(((xt1 - xr1) * dx1) * ((xt1 - xr1) * dx1) + ((yt1 - yr1) * dy1) * ((yt1 - yr1) * dy1));// ͨ����·����
		N_Sample = (int)floor(Length_L / acc);		// ��·�ܵļ������
		spacX = (xr1 - xt1) / N_Sample;		// �������������ľ��Ȳ�ֵ
		spacY = (yr1 - yt1) / N_Sample;		// ��������������γ�Ȳ�ֵ
		Fvalue = new double*[1];
		Fvalue[0] = new double[N_Sample];
		if (Fvalue == NULL)
		{
			cout<<"�ڴ治�㣡";
		}

		////////////////////////// ITU-rp526��ǿԤ��ģ���� //////////////////////////////
		ITUrp RaosheITU;
		for (j = 0; j < N_Sample; j++)
		{
			xi = xt1 + floor(spacX * (j + 1));
			yj = yt1 + floor(spacY * (j + 1));
			lined = sqrt(((xt1 - xi) * dx1) * ((xt1 - xi) * dx1) + ((yt1 - yj) * dy1) * ((yt1 - yj) * dy1));
			if ((xi == floor(xt1)) && (yj == floor(yt1)))
			{
				Fvalue[0][j] = 0.0;
			}
			//else if (lined > distance * 1000)// ��������Ϊ����ͨ��
			//{
			//	Fvalue[0][j] = 1000;	// ����ͨ�ŷ�Χ����Ϊ��ļ�����Ϊ˥��1000dB
			//}
			/*else if (lined < 1000)
			{
			Fvalue[0][j] = 32.45 + 20 * log10(lined / 1000) + 26 * log10(f);
			}*/
			else if (lined < 300)
			{
				Fvalue[0][j] = 32.45 + 20 * log10(lined / 1000) + 26 * log10(f);
			}
			else
			{
				Fvalue[0][j] = RaosheITU.ITUraoshe(xt1, yt1, xi, yj, GRID, row, column, acc, f, ht, hr, y3, y4, aa, e, re, dx1, dy1);
			}
		}
	}




	for (i = 0; i<row; i++)
	{
		delete[]GRID[i];
	}
	delete[]GRID;
	return Fvalue;
}


double **DecayDiffraction(double Setparameter[],double Receparameter[],Data **Eledata, int dimension[],double f,Iparameter IPara)
{
	double xt,yt,Gt, ht,xr,yr,Gr,hr, acc, aa, distance;
	int row, column, e, i, j;
	xt=Setparameter[0]; // �������߾���
	yt=Setparameter[1]; // ��������γ��
	Gt=Setparameter[2]; // �������������������
	ht=Setparameter[3]; // �������߸߶�

	xr=Receparameter[0]; // �������߾���
	yr=Receparameter[1]; // ��������γ��
	Gr=Receparameter[2]; // �������������������
	hr=Receparameter[3]; // �������߸߶�

	row = dimension[0];    // ������Ϣ���ݵ�����
	column = dimension[1]; // ������Ϣ���ݵ�����
	acc = IPara.pSize;	   // Ԥ�������
	distance = IPara.Elength;// ��ЧԤ�ⷶΧ
	if(IPara.aMark == 0)// ��ˮ
	{
		aa=0.001;		// ��Ч�絼��
		e=80;			// ��Խ�糣��
	}
	if(IPara.aMark == 1)// ��ˮ
	{
		aa=4.0;			// ��Ч�絼��
		e=80;			// ��Խ�糣��
	}
	if(IPara.aMark == 2)// ʪ��
	{
		aa=0.008;		// ��Ч�絼��
		e=10;			// ��Խ�糣��
	}
	if(IPara.aMark == 3)// �ɵ�
	{
		aa=0.001;		// ��Ч�絼��
		e=2;			// ��Խ�糣��
	}
	double **Fvalue; // ���巵��ֵ��άָ��
	double **GRID;   // ����߳̾���ָ��
	double dx,dy,dx1,dy1,xt1,yt1,xr1,yr1;
	int y3=1,y4=0;		// ���������͵������
	double re=8500000;  // �����Ч�뾶
	dy = fabs(Eledata[0][0].ycor-Eledata[1][1].ycor); // ���������������γ�Ȳ�ֵ
	dx = fabs(Eledata[0][0].xcor-Eledata[1][1].xcor); // ��������������ľ��Ȳ�ֵ
	// ======================== ��ͼ�����ڲ������ˮƽ���ʹ�ֱ���m ============================ //
	//dx1 = dx * 6371000 * cos(Eledata[0][1].ycor/180*3.1415926) * 3.1415926/180.0; // ��ͼˮƽ������m
	//dy1 = dy * 6371000 * 3.1415926/180.0;                                         // ��ͼ��ֱ������m
    dx1 = 73.5624;			  // ��ͼˮƽ������m
	dy1 = 92.6043;			  // ��ͼ��ֱ������m
	// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = //
	xt1 = floor(fabs(xt - Eledata[0][0].xcor)/dx);    // ������ڵ�ͼ�����ϵ�λ��
	yt1 = floor(fabs(yt - Eledata[0][0].ycor)/dy);
	xr1 = floor(fabs(xr - Eledata[0][0].xcor)/dx);    // ���յ��ڵ�ͼ�����ϵ�λ��
	yr1 = floor(fabs(yr - Eledata[0][0].ycor)/dy);	
	GRID = new double*[row];
	for(i=0;i < row;i++)
	{
		GRID[i] = new double[column];
	}
	for(i = 0; i < row; i++)
	{
		for(j = 0;j < column; j++)
		{
			GRID[i][j] = Eledata[i][j].eleva;
		}
	}
	

		// ��̬�����ڴ�
		Fvalue = new double*[row];
		for(i=0;i < row;i++)
		{
			Fvalue[i] = new double[column]; 
		}
		//////////////////////// ITU-rp526��ǿԤ��ģ���� ////////////////////////////////////////////
		ITUrp RaosheITU; 
		double lined;
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < column; j++)
			{
				lined = sqrt(((xt1 - i) * dx1) * ((xt1 - i) * dx1) + ((yt1 - j) * dy1) * ((yt1 - j) * dy1)); //����վ����ͼĳ��ľ���
				if (lined == 0)
				{
					Fvalue[i][j] = 0.0;
				}
				else if (lined > distance*1000)// ��������Ϊ����ͨ��
				{
					Fvalue[i][j] = 1000;	// ����ͨ�ŷ�Χ����Ϊ��ļ�����Ϊ˥��1000dB
				}
				else if (lined < 1000)
				{
					Fvalue[i][j] = 32.45 + 20 * log10(lined / 1000) + 26 * log10(f); //���ɿռ䴫��
				}
				else
				{
					Fvalue[i][j] = RaosheITU.ITUraoshe(xt1,yt1,i,j,GRID,row,column,acc,f,ht,0,y3,y4,aa, e, re, dx1, dy1);
				}
			}
		}

	
	
	for(i=0; i<row;i++)
	{
		delete[]GRID[i];
	}
	delete []GRID;
	return Fvalue;
}

double cauclateGain(double xt,double yt,double xr,double yr,double Azimuth,double diameter,double fre)
{
	double GGr;
	double arfa,arfaA,phai,phai1,phair;
	double lamde = 3.0/fre;
	double Gain0 = 10*log10(0.6*(PI*diameter/lamde)*(PI*diameter/lamde));

	arfa = atan(sin(xr-xt)/(cos(yt)*tan(yr)-sin(yt)*cos(xr-xt)));
	if(tan(arfa > 0))
	{
		if(yt < yr)
		{
			arfaA = arfa;
		}
		else
		{
			arfaA = PI+arfa;		
		}
	}
	else
	{
		if(yt < yr)
		{
			arfaA = 2*PI-fabs(arfa);
		}
		else
		{
			arfaA = PI-fabs(arfa);		
		}
	}

	phai = fabs(Azimuth - arfaA)/PI*180;
	phai1 = 20*lamde*sqrt(Gain0-2-15*log10(diameter/lamde))/diameter;
	phair = 15.85*pow((diameter/lamde),-0.6);
	if(phai<phai1)
	{
		GGr = Gain0-2.5*(phai*diameter/lamde)*(phai*diameter/lamde)/1000;
	}
	else if(phai<100*lamde/diameter)
	{
		GGr = 2+15*log10(diameter/lamde);
	}
	else if(phai<48)
	{
		GGr = 52-10*log10(diameter/lamde)-25*log10(phai);
	}
	else
	{
		GGr = 10 - 10*log(diameter/lamde);
	}
	return GGr;
}

double *DecayMicrowave(double Setparameter[],double Receparameter[],Data **Eledata, int dimension[],double f,ParWave Para)
{
	double xt,yt,Gt,ht,Azimutht,xr,yr,Gr,hr,Azimuthr, acc, distance, zetaN, PW, Diameter;
	int row, column, i, j;
	xt=Setparameter[0]; // �������߾���
	yt=Setparameter[1]; // ��������γ��
	Gt=Setparameter[2]; // �������������������
	ht=Setparameter[3]; // �������߸߶�
	Azimutht = Setparameter[4]/180*PI; // �������߷�λ��

	xr=Receparameter[0]; // �������߾���
	yr=Receparameter[1]; // ��������γ��
	Gr=Receparameter[2]; // �������������������
	hr=Receparameter[3]; // �������߸߶�
	Azimuthr = Receparameter[4]/180*PI;// �������߷�λ��

	row = dimension[0];    // ������Ϣ���ݵ�����
	column = dimension[1]; // ������Ϣ���ݵ�����
	acc = Para.pSize*1000; // Ԥ������� ��λm
	distance = Para.Elength;// ��ЧԤ�ⷶΧ
	zetaN = Para.deterN;	// ����������ָ���ݶ�
	PW = Para.pw;			// ����·ݵ�ʱ��ٷֱ�
	Diameter = Para.diameter;// ����ֱ��
	double fre = f/1000;
	double *Fvalue; // ���巵��ֵһάָ��
	double **GRID;   // ����߳̾���ָ��
	double dx,dy,dx1,dy1,xt1,yt1,xr1,yr1,Gaint;
	double re=8500000;  // �����Ч�뾶
	double fai = (yt+yr)/2.0;
	dy = fabs(Eledata[0][0].ycor-Eledata[1][1].ycor); // ���������������γ�Ȳ�ֵ
	dx = fabs(Eledata[0][0].xcor-Eledata[1][1].xcor); // ��������������ľ��Ȳ�ֵ
	// ======================== ��ͼ�����ڲ������ˮƽ���ʹ�ֱ���m ============================ //
	//dx1 = dx * 6371000 * cos(Eledata[0][1].ycor/180*3.1415926) * 3.1415926/180.0; // ��ͼˮƽ������m
	//dy1 = dy * 6371000 * 3.1415926/180.0;                                         // ��ͼ��ֱ������m
    dx1 = 73.5624;			  // ��ͼˮƽ������m
	dy1 = 92.6043;			  // ��ͼ��ֱ������m
	// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = //
	xt1 = floor(fabs(xt - Eledata[0][0].xcor)/dx);    // ������ڵ�ͼ�����ϵ�λ��
	yt1 = floor(fabs(yt - Eledata[0][0].ycor)/dy);
	xr1 = floor(fabs(xr - Eledata[0][0].xcor)/dx);    // ���յ��ڵ�ͼ�����ϵ�λ��
	yr1 = floor(fabs(yr - Eledata[0][0].ycor)/dy);	
	GRID = new double*[row];
	for(i=0;i < row;i++)
	{
		GRID[i] = new double[column];
	}
	for(i = 0; i < row; i++)
	{
		for(j = 0;j < column; j++)
		{
			GRID[i][j] = Eledata[i][j].eleva;
		}
	}
	// = = = = = = = = = = = = = = = = = = = = ΢����·˥��Ԥ��ģ�� = = = = = = = = = = = = = = = = = = //		
	double lamde = 3.0/fre;
	//Gaint = 10*log10(0.6*(PI*Diameter/lamde)*(PI*Diameter/lamde));
	Gaint = cauclateGain(xt,yt,xr,yr,Azimutht,Diameter,fre);
	double Length_L,lined,spacX,spacY,xi,yj;
	int N_Sample;
	Length_L = sqrt(((xt1 - xr1) * dx1) * ((xt1 - xr1) * dx1) + ((yt1 - yr1) * dy1) * ((yt1 - yr1) * dy1));// ͨ����·����
	N_Sample = (int)floor(Length_L/acc);		// ��·�ܵļ������
	spacX = ((xt1 - xr1) * dx1)/N_Sample;	// �������������ľ��Ȳ�ֵ
	spacY = ((yt1 - yr1) * dy1)/N_Sample;	// ��������������γ�Ȳ�ֵ
	Fvalue = new double[N_Sample]; 
	attenuation_ra raITU;
	for(j = 0;j < N_Sample; j++)
	{
		xi = xt1 + floor(spacX * (j + 1) / dx1);
		yj = yt1 + floor(spacY * (j + 1) / dy1);
		lined = sqrt(((xt1 - xi) * dx1) * ((xt1 - xi) * dx1) + ((yt1 - yj) * dy1) * ((yt1 - yj) * dy1));
		if ((xi == floor(xt1)) && (yj == floor(yt1)))
		{
			Fvalue[j] = 0.0;
		}
		else if (lined > distance*1000)// ��������Ϊ����ͨ��
		{
			Fvalue[j] = 1000;	// ����ͨ�ŷ�Χ����Ϊ��ļ�����Ϊ˥��1000dB
		}
		else if (lined < 1000)
		{
			Fvalue[j] = 92.5 + 20 * log10(lined / 1000) + 20 * log10(fre);
		}
		else
		{
			Fvalue[j] = raITU.Radio_C(GRID,zetaN,xt1,yt1,xi,yj,row,column,fre,Gaint,0,Para.pSize,dx1,dy1,PW,fai,ht);
		}
	}		

	for(i=0; i<row;i++)
	{
		delete[]GRID[i];
	}
	delete []GRID;
	return Fvalue;
}

double PointDiffraction(double Setparameter[],double Receparameter[],Data **Eledata, int dimension[],double f,Iparameter IPara)
{
	double xt,yt,Gt, ht,xr,yr,Gr,hr, acc, aa, distance;
	int row, column, e, i, j;
	xt=Setparameter[0]; // �������߾���
	yt=Setparameter[1]; // ��������γ��
	Gt=Setparameter[2]; // �������������������
	ht=Setparameter[3]; // �������߸߶�

	xr=Receparameter[0]; // �������߾���
	yr=Receparameter[1]; // ��������γ��
	Gr=Receparameter[2]; // �������������������
	hr=Receparameter[3]; // �������߸߶�

	row = dimension[0];    // ������Ϣ���ݵ�����
	column = dimension[1]; // ������Ϣ���ݵ�����
	acc = IPara.pSize;	   // Ԥ�������
	distance = IPara.Elength;// ��ЧԤ�ⷶΧ
	if(IPara.aMark == 0)// ��ˮ
	{
		aa=0.001;		// ��Ч�絼��
		e=80;			// ��Խ�糣��
	}
	if(IPara.aMark == 1)// ��ˮ
	{
		aa=4.0;			// ��Ч�絼��
		e=80;			// ��Խ�糣��
	}
	if(IPara.aMark == 2)// ʪ��
	{
		aa=0.008;		// ��Ч�絼��
		e=10;			// ��Խ�糣��
	}
	if(IPara.aMark == 3)// �ɵ�
	{
		aa=0.001;		// ��Ч�絼��
		e=2;			// ��Խ�糣��
	}
	double ret_Value; // ���巵��ֵ
	double **GRID;    // ����߳̾���ָ��
	double dx,dy,dx1,dy1,xt1,yt1,xr1,yr1;
	int y3=1,y4=0;		// ���������͵������
	double re=8500000;  // �����Ч�뾶
	dy = fabs(Eledata[0][0].ycor-Eledata[1][1].ycor); // ���������������γ�Ȳ�ֵ
	dx = fabs(Eledata[0][0].xcor-Eledata[1][1].xcor); // ��������������ľ��Ȳ�ֵ
	// ======================== ��ͼ�����ڲ������ˮƽ���ʹ�ֱ���m ============================ //
	//dx1 = dx * 6371000 * cos(Eledata[0][1].ycor/180*3.1415926) * 3.1415926/180.0; // ��ͼˮƽ������m
	//dy1 = dy * 6371000 * 3.1415926/180.0;                                         // ��ͼ��ֱ������m
    dx1 = 73.5624;			  // ��ͼˮƽ������m
	dy1 = 92.6043;			  // ��ͼ��ֱ������m
	// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = //
	xt1 = floor(fabs(xt - Eledata[0][0].xcor)/dx);    // ������ڵ�ͼ�����ϵ�λ��
	yt1 = floor(fabs(yt - Eledata[0][0].ycor)/dy);
	xr1 = floor(fabs(xr - Eledata[0][0].xcor)/dx);    // ���յ��ڵ�ͼ�����ϵ�λ��
	yr1 = floor(fabs(yr - Eledata[0][0].ycor)/dy);	
	GRID = new double*[row];
	for(i=0;i < row;i++)
	{
		GRID[i] = new double[column];
	}
	for(i = 0; i < row; i++)
	{
		for(j = 0;j < column; j++)
		{
			GRID[i][j] = Eledata[i][j].eleva;
		}
	}
	// = = = = = = = = = = = = = = = = = = = = ITU-rp526ģ��˥��Ԥ�� = = = = = = = = = = = = = = = = = = //
	double lined;
	lined = sqrt(((xt1 - xr1) * dx1) * ((xt1 - xr1) * dx1) + ((yt1 - yr1) * dy1) * ((yt1 - yr1) * dy1));// ͨ����·����
	ITUrp RaosheITU;// ITU-rp526��ǿԤ��ģ���� 
	if (lined == 0)
	{
		ret_Value = 0.0;
	}
	else if (lined > distance*1000)// ��������Ϊ����ͨ��
	{
		ret_Value = 1000;	// ����ͨ�ŷ�Χ����Ϊ��ļ�����Ϊ˥��1000dB
	}
	else if (lined < 1000)
	{
		ret_Value = 32.45 + 20 * log10(lined / 1000) + 26 * log10(f);
	}
	else
	{
		ret_Value = RaosheITU.ITUraoshe(xt1,yt1,xr1,yr1,GRID,row,column,acc,f,ht,hr,y3,y4,aa, e, re, dx1, dy1);
	}
		
	// = = = = = = = = = = = = = = = = = = = = �ͷ��ڴ� = = = = = = = = = = = = = = = = = = = = = = = = = //	
	for(i=0; i<row;i++)
	{
		delete[]GRID[i];
	}
	delete []GRID;

	return ret_Value;
}


double PointMicrowave(double Setparameter[],double Receparameter[],Data **Eledata, int dimension[],double f,ParWave Para)
{
	double xt,yt,Gt,ht,Azimutht,xr,yr,Gr,hr,Azimuthr, acc, distance, zetaN, PW, Diameter;
	int row, column, i, j;
	xt=Setparameter[0]; // �������߾���
	yt=Setparameter[1]; // ��������γ��
	Gt=Setparameter[2]; // �������������������
	ht=Setparameter[3]; // �������߸߶�
	Azimutht = Setparameter[4]/180*PI; // �������߷�λ��

	xr=Receparameter[0]; // �������߾���
	yr=Receparameter[1]; // ��������γ��
	Gr=Receparameter[2]; // �������������������
	hr=Receparameter[3]; // �������߸߶�
	Azimuthr = Receparameter[4]/180*PI;// �������߷�λ��

	row = dimension[0];    // ������Ϣ���ݵ�����
	column = dimension[1]; // ������Ϣ���ݵ�����
	acc = Para.pSize*1000; // Ԥ������� ��λm
	distance = Para.Elength;// ��ЧԤ�ⷶΧ
	zetaN = Para.deterN;	// ����������ָ���ݶ�
	PW = Para.pw;			// ����·ݵ�ʱ��ٷֱ�
	Diameter = Para.diameter;// ����ֱ��
	double fre = f/1000;
	double ret_Value; // ���巵��ֵһάָ��
	double **GRID;   // ����߳̾���ָ��
	double dx,dy,dx1,dy1,xt1,yt1,xr1,yr1,Gaint,Gainr;
	double re=8500000;  // �����Ч�뾶
	double height1 = ht; 
	double fai = (yt+yr) / 2.0;
	dy = fabs(Eledata[0][0].ycor-Eledata[1][1].ycor); // ���������������γ�Ȳ�ֵ
	dx = fabs(Eledata[0][0].xcor-Eledata[1][1].xcor); // ��������������ľ��Ȳ�ֵ
	// ======================== ��ͼ�����ڲ������ˮƽ���ʹ�ֱ���m ============================ //
	//dx1 = dx * 6371000 * cos(Eledata[0][1].ycor/180*3.1415926) * 3.1415926/180.0; // ��ͼˮƽ������m
	//dy1 = dy * 6371000 * 3.1415926/180.0;                                         // ��ͼ��ֱ������m
    dx1 = 73.5624;			  // ��ͼˮƽ������m
	dy1 = 92.6043;			  // ��ͼ��ֱ������m
	// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = //
	xt1 = floor(fabs(xt - Eledata[0][0].xcor)/dx);    // ������ڵ�ͼ�����ϵ�λ��
	yt1 = floor(fabs(yt - Eledata[0][0].ycor)/dy);
	xr1 = floor(fabs(xr - Eledata[0][0].xcor)/dx);    // ���յ��ڵ�ͼ�����ϵ�λ��
	yr1 = floor(fabs(yr - Eledata[0][0].ycor)/dy);	
	GRID = new double*[row];
	for(i=0;i < row;i++)
	{
		GRID[i] = new double[column];
	}
	for(i = 0; i < row; i++)
	{
		for(j = 0;j < column; j++)
		{
			GRID[i][j] = Eledata[i][j].eleva;
		}
	}
	// = = = = = = = = = = = = = = = = = = = = ΢����·˥��Ԥ��ģ�� = = = = = = = = = = = = = = = = = = //		
	double lamde = 3.0/fre;
	//Gaint = 10*log10(0.6*(PI*Diameter/lamde)*(PI*Diameter/lamde));
	Gaint = cauclateGain(xt,yt,xr,yr,Azimutht,Diameter,fre);
	Gainr = cauclateGain(xr,yr,xt,yt,Azimuthr,Diameter,fre);
	double lined;
	lined = sqrt(((xt1 - xr1) * dx1) * ((xt1 - xr1) * dx1) + ((yt1 - yr1) * dy1) * ((yt1 - yr1) * dy1));// ͨ����·����
	
	attenuation_ra raITU;
	if (lined == 0)
	{
		ret_Value = 0.0;
	}
	else if (lined > distance*1000)// ��������Ϊ����ͨ��
	{
		ret_Value = 1000;	// ����ͨ�ŷ�Χ����Ϊ��ļ�����Ϊ˥��1000dB
	}
	else if (lined < 1000)
	{
		ret_Value = 32.45 + 20 * log10(lined / 1000) + 26 * log10(f);
	}
	else
	{
		ret_Value = raITU.Radio_C(GRID,zetaN,xt1,yt1,xr1,yr1,row,column,fre,Gaint,Gainr,Para.pSize,dx1,dy1,PW,fai,height1);
	}
	

	for(i=0; i<row;i++)
	{
		delete[]GRID[i];
	}
	delete []GRID;
	return ret_Value;
}
