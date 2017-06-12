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
	double xcor;    // 经度
	double ycor;    // 纬度
	double eleva;   // 高程值
};
struct Iparameter
{
	double pSize;   // 预测颗粒度
	double Elength; // 有效预测范围
	int aMark;      // 电导率和相对介电常数标志
};
typedef double**(*Calculate)(double[],double[],Data **,int[],double,Iparameter,double,double,int);

int main(int argc, char* argv[])
{
		//========================= 地图文件的读取 =================================//
    double DEM_x0, DEM_y0, cells, tmp;
	char ncols[20],nrows[20],xcorner[20],ycorner[20],csize[20],NODATA[20];
	int row,column, value, i, j;
	FILE *fp;
	// 打开目标地图文件
	fp=fopen("N44E125.asc","r");	
	if(fp==NULL)
	{
		cout<< "can't open the GRD datafile!" <<endl;
		return 0;
	}
	// 读取文件头信息
	fscanf(fp,"%s %d %s %d %s %lf %s %lf %s %lf %s %lf",
            &ncols, &column, &nrows, &row, &xcorner, &DEM_x0,
            &ycorner, &DEM_y0, &csize, &cells, &NODATA,&value);
	printf("输出地图左下角的经、纬度值：\n");
	printf("%f  ",DEM_x0);
	printf("%f\n",DEM_y0);
	printf("输出高程值的行列数：\n");
	printf("%d  ",row);
	printf("%d\n",column);
	// 读取文件高程数据实体
	double **DEM;
	DEM = new double*[row];
	for (int i = 0; i < column; i++)
	{
		DEM[i] = new double[column];
	}
	if (DEM == NULL)
	{
		cout << "内存不足！";
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

	//============================ 初始化信息 =====================================//
	double xt = 125.625416666633, yt = 44.5137499999666, Gt = 15, h_g1 = 15;      // 发射端初始参数
	double xr = 125.8070833333, yr = 44.3554166666333, Gr = 15, h_g2 = 15;     // 接收端初始参数
	double f = 60; // 电磁波频率，单位MHz
	int mark = 1;		// 标志位：0表示区域覆盖预测；1表示传输路径损耗预测
	double dx = cells, dy = cells;	// 地图水平、垂直两点的间距（度）
	Data **EleData;
	EleData = new Data*[row];  // 动态申请内存
	for (i = 0; i < row; i++)
	{
		EleData[i] = new Data[column];
	}
	for (i = 0; i < row; i++)
	{
		for (j = 0; j < column; j++)
		{
			EleData[i][j].eleva = DEM[i][j]; //高程值
			EleData[i][j].xcor = DEM_x0 + dx * j;    // 经度
			EleData[i][j].ycor = DEM_y0 + dy * (column - i);//纬度
		}
	}
	double dx1 = 73.5624;			  // 地图水平两点间距m
	double dy1 = 92.6043;			  // 地图垂直两点间距m 1200*1200的地图
	//double dx1 = 882.7488;			 //100*100的地图
	//double dy1 = 1111.2516;
	//double dx1 = 276.0889;			 //32*32的地图,0.1个经度
	//double dy1 = 347.5555;
	/*double dx1 = 5;
	double dy1 = 5;*/

	
	int xt1, yt1, xr1, yr1;
	xt1 = (int)floor(fabs(xt - EleData[0][0].xcor) / dx);    // 发射点在地图矩阵上的位置
	yt1 = (int)floor(fabs(yt - EleData[0][0].ycor) / dy);
	xr1 = (int)floor(fabs(xr - EleData[0][0].xcor) / dx);    // 接收点在地图矩阵上的位置
	yr1 = (int)floor(fabs(yr - EleData[0][0].ycor) / dy);
	//xt1 = 48;   //侧面遮挡测试
	//yt1 = 1;
	//xr1 = 48;
	//yr1 = 95;

	double Setparameter[4],Receparameter[4];
	Setparameter[0]=xt1;			// 发射点经度
	Setparameter[1]=yt1;			// 发射点纬度
	Setparameter[2]=Gt;			// 发射天线增益
	Setparameter[3]=h_g1;		// 发射天线高度
	Receparameter[0]=xr1;		// 接收点经度
	Receparameter[1]=yr1;		// 接收点纬度
	Receparameter[2]=Gr;		// 接收天线增益
	Receparameter[3]=h_g2;		// 接收天线高度
	
	int dimension[2];
	dimension[0]=row;			// 地理高程矩阵的行数
	dimension[1]=column;		// 地理高程矩阵的列数
	Iparameter IPara;
	IPara.aMark = 2;			// 电导率和相对介电常数标志
	IPara.Elength = 30;			// 有效通信范围，单位km
	IPara.pSize = 150;			// 计算精度
	
	//for (i = 0; i<row; i++)
	//{
	//	delete[]DEM[i];
	//}
	//delete[]DEM;
	 //========================== 与输出相关的计算 =================================//

	double LengLine=sqrt((xt1-xr1)*dx1*(xt1-xr1)*dx1+(yt1-yr1)*dy1*(yt1-yr1)*dy1);
	int N_Sample = (int)floor(LengLine/IPara.pSize);			   // 链路总的计算点数
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


	// 调取动态链接库计算出衰减矩阵或衰减数组
	HINSTANCE hDLL;
	Calculate CountFun;
	hDLL = LoadLibrary(".\\ForecastDecay.dll");
	if(hDLL != NULL)
	{
		CountFun = (Calculate)GetProcAddress(hDLL,"DecayDiffraction");
		if(CountFun != NULL)
		{
			Fvalue = CountFun(Setparameter,Receparameter,EleData,dimension,f,IPara,dx1,dy1,mark);
			// 输出预测结果
			if(mark ==0)
			{
				cout << "区域覆盖预测结果：" << endl;
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
				printf("传输路径损耗预测结果：\n");
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