// ElecCall.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <stdio.h>
#include <windows.h>
#include <math.h>
#include <fstream>

using namespace std;

struct Data
{
	double xcor;    // ����
	double ycor;    // γ��
	double eleva;   // �߳�ֵ
};
struct Iparameter
{
	double pSize;   // Ԥ�������
	double Elength; // ��ЧԤ�ⷶΧ
	int aMark;      // �絼�ʺ���Խ�糣����־
};
typedef double**(*Calculate)(double[],double[],Data **,int[],double,Iparameter,double,double,int);

int main(int argc, char* argv[])
{
		//========================= ��ͼ�ļ��Ķ�ȡ =================================//
    double DEM_x0, DEM_y0, cells, tmp;
	char ncols[20],nrows[20],xcorner[20],ycorner[20],csize[20],NODATA[20];
	int row,column, value, i, j;
	FILE *fp;
	// ��Ŀ���ͼ�ļ�
	fp=fopen("N44E125.asc","r");	
	if(fp==NULL)
	{
		cout<< "can't open the GRD datafile!" <<endl;
		return 0;
	}
	// ��ȡ�ļ�ͷ��Ϣ
	fscanf(fp,"%s %d %s %d %s %lf %s %lf %s %lf %s %lf",
            &ncols, &column, &nrows, &row, &xcorner, &DEM_x0,
            &ycorner, &DEM_y0, &csize, &cells, &NODATA,&value);
	printf("�����ͼ���½ǵľ���γ��ֵ��\n");
	printf("%f  ",DEM_x0);
	printf("%f\n",DEM_y0);
	printf("����߳�ֵ����������\n");
	printf("%d  ",row);
	printf("%d\n",column);
	// ��ȡ�ļ��߳�����ʵ��
	double **DEM;
	DEM = new double*[row];
	for (int i = 0; i < column; i++)
	{
		DEM[i] = new double[column];
	}
	if (DEM == NULL)
	{
		cout << "�ڴ治�㣡";
		return 0;
	}
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < column; j++)
		{
			fscanf(fp, "%lf", &tmp);
			if (tmp == -9999)
				DEM[i][j] = -9999;
			else
				DEM[i][j] = tmp;
		}

	}
	fclose(fp);

	//============================ ��ʼ����Ϣ =====================================//
	double xt = 125.625416666633, yt = 44.5137499999666, Gt = 15, h_g1 = 15;      // ����˳�ʼ����
	double xr = 125.8070833333, yr = 44.3554166666333, Gr = 15, h_g2 = 15;     // ���ն˳�ʼ����
	double f = 60; // ��Ų�Ƶ�ʣ���λMHz
	int mark = 1;		// ��־λ��0��ʾ���򸲸�Ԥ�⣻1��ʾ����·�����Ԥ��
	double dx = cells, dy = cells;	// ��ͼˮƽ����ֱ����ļ�ࣨ�ȣ�
	Data **EleData;
	EleData = new Data*[row];  // ��̬�����ڴ�
	for (i = 0; i < row; i++)
	{
		EleData[i] = new Data[column];
	}
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < column; j++)
		{
			EleData[i][j].eleva = DEM[i][j]; //�߳�ֵ
			EleData[i][j].xcor = DEM_x0 + dx * j;    // ����
			EleData[i][j].ycor = DEM_y0 + dy * (column - i);//γ��
		}
	}
	double dx1 = 73.5624;			  // ��ͼˮƽ������m
	double dy1 = 92.6043;			  // ��ͼ��ֱ������m 1200*1200�ĵ�ͼ
	//double dx1 = 882.7488;			 //100*100�ĵ�ͼ
	//double dy1 = 1111.2516;
	//double dx1 = 276.0889;			 //32*32�ĵ�ͼ,0.1������
	//double dy1 = 347.5555;
	/*double dx1 = 5;
	double dy1 = 5;*/

	
	int xt1, yt1, xr1, yr1;
	xt1 = (int)floor(fabs(xt - EleData[0][0].xcor) / dx);    // ������ڵ�ͼ�����ϵ�λ��
	yt1 = (int)floor(fabs(yt - EleData[0][0].ycor) / dy);
	xr1 = (int)floor(fabs(xr - EleData[0][0].xcor) / dx);    // ���յ��ڵ�ͼ�����ϵ�λ��
	yr1 = (int)floor(fabs(yr - EleData[0][0].ycor) / dy);
	//xt1 = 48;   //�����ڵ�����
	//yt1 = 1;
	//xr1 = 48;
	//yr1 = 95;

	double Setparameter[4],Receparameter[4];
	Setparameter[0]=xt1;			// ����㾭��
	Setparameter[1]=yt1;			// �����γ��
	Setparameter[2]=Gt;			// ������������
	Setparameter[3]=h_g1;		// �������߸߶�
	Receparameter[0]=xr1;		// ���յ㾭��
	Receparameter[1]=yr1;		// ���յ�γ��
	Receparameter[2]=Gr;		// ������������
	Receparameter[3]=h_g2;		// �������߸߶�
	
	int dimension[2];
	dimension[0]=row;			// ����߳̾��������
	dimension[1]=column;		// ����߳̾��������
	Iparameter IPara;
	IPara.aMark = 2;			// �絼�ʺ���Խ�糣����־
	IPara.Elength = 30;			// ��Чͨ�ŷ�Χ����λkm
	IPara.pSize = 150;			// ���㾫��
	
	//for (i = 0; i<row; i++)
	//{
	//	delete[]DEM[i];
	//}
	//delete[]DEM;
	 //========================== �������صļ��� =================================//

	double LengLine=sqrt((xt1-xr1)*dx1*(xt1-xr1)*dx1+(yt1-yr1)*dy1*(yt1-yr1)*dy1);
	int N_Sample = (int)floor(LengLine/IPara.pSize);			   // ��·�ܵļ������
	double **Fvalue;
	if(mark ==0)
	{
		Fvalue = new double*[row];
		for(i=0;i < row;i++)
		{
			Fvalue[i] = new double[column]; 
		}
	}
	else
	{
		Fvalue = new double*[1];
		Fvalue[i] = new double[N_Sample]; 
	}


	ofstream f1("..\\t1.txt");
	if (!f1)return 0;


	// ��ȡ��̬���ӿ�����˥�������˥������
	HINSTANCE hDLL;
	Calculate CountFun;
	hDLL = LoadLibrary(".\\ForecastDecay.dll");
	if(hDLL != NULL)
	{
		CountFun = (Calculate)GetProcAddress(hDLL,"DecayDiffraction");
		if(CountFun != NULL)
		{
			Fvalue = CountFun(Setparameter,Receparameter,EleData,dimension,f,IPara,dx1,dy1,mark);
			// ���Ԥ����
			if(mark ==0)
			{
				cout << "���򸲸�Ԥ������" << endl;
				for (i = 0; i < row; i++)
				{
					for (j = 0; j < column; j++)
					{
						if (j != column - 1)
						{
							//printf("%f ", Fvalue[i][j]);
							f1 << Fvalue[i][j] << " ";
						}
						else
						{
							//printf("%f\n", Fvalue[i][j]);
							f1 << Fvalue[i][j] << endl;
						}


					}
				}
				f1.close();
			}
			else
			{
				printf("����·�����Ԥ������\n");
				for(j = 0;j < N_Sample; j++)
				{
					f1 << Fvalue[0][j] << " ";
				}
			}		
		}
		
		FreeLibrary(hDLL);
	}
	


	system("pause");
	return 0;
}