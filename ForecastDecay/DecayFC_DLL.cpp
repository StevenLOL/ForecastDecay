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
	xt1 = Setparameter[0]; // 发射天线经度
	yt1 = Setparameter[1]; // 发射天线纬度
	Gt = Setparameter[2]; // 发射天线最大天线增益
	ht = Setparameter[3]; // 发射天线高度

	xr1 = Receparameter[0]; // 接收天线经度
	yr1 = Receparameter[1]; // 接收天线纬度
	Gr = Receparameter[2]; // 接收天线最大天线增益
	hr = Receparameter[3]; // 接收天线高度

	row = dimension[0];    // 地理信息数据的行数
	column = dimension[1]; // 地理信息数据的列数
	acc = IPara.pSize;	   // 预测颗粒度
	distance = IPara.Elength;// 有效预测范围
	if (IPara.aMark == 0)// 淡水
	{
		aa = 0.001;		// 有效电导率
		e = 80;			// 相对介电常数
	}
	if (IPara.aMark == 1)// 海水
	{
		aa = 4.0;			// 有效电导率
		e = 80;			// 相对介电常数
	}
	if (IPara.aMark == 2)// 湿地
	{
		aa = 0.008;		// 有效电导率
		e = 10;			// 相对介电常数
	}
	if (IPara.aMark == 3)// 干地
	{
		aa = 0.001;		// 有效电导率
		e = 2;			// 相对介电常数
	}
	double **Fvalue; // 定义返回值二维指针
	double **GRID;   // 定义高程矩阵指针
	int y3 = 1, y4 = 0;		// 极化参数和地面参数
	double re = 8500000;  // 地球等效半径
	// ======================== 地图中相邻采样点的水平间距和垂直间距m ============================ //

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

 //mark=0  区域场强预测
	if (mark == 0)
	{
		Fvalue = new double*[row];
			for (i = 0; i < row; i++)
			{
				Fvalue[i] = new double[column];
			}
			//////////////////////// ITU-rp526场强预测模型类 ////////////////////////////////////////////


			ITUrp RaosheITU;
			double lined;
			for (i = 0; i < row; i++)
			{
				for (j = 0; j < column; j++)
				{
					lined = sqrt(((xt1 - i) * dx1) * ((xt1 - i) * dx1) + ((yt1 - j) * dy1) * ((yt1 - j) * dy1));//发射站到地图某点的距离
					if (lined == 0)
					{
						Fvalue[i][j] = 0.0;
					}
					//else if (lined > distance * 1000)// 超出，认为不可通信
					//{
					//	Fvalue[i][j] = 1000;	// 超出通信范围，认为损耗极大，设为衰减1000db
					//}
					else if (lined < 1000)
					{
						Fvalue[i][j] = 32.45 + 20 * log10(lined / 1000) + 26 * log10(f); //自由空间传播
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
		//mark=1 链路衰减预测
		double Length_L, lined, spacX, spacY, xi, yj;
		int N_Sample;
		Length_L = sqrt(((xt1 - xr1) * dx1) * ((xt1 - xr1) * dx1) + ((yt1 - yr1) * dy1) * ((yt1 - yr1) * dy1));// 通信链路长度
		N_Sample = (int)floor(Length_L / acc);		// 链路总的计算点数
		spacX = (xr1 - xt1) / N_Sample;		// 相邻两个计算点的经度差值
		spacY = (yr1 - yt1) / N_Sample;		// 相邻两个计算点的纬度差值
		Fvalue = new double*[1];
		Fvalue[0] = new double[N_Sample];
		if (Fvalue == NULL)
		{
			cout<<"内存不足！";
		}

		////////////////////////// ITU-rp526场强预测模型类 //////////////////////////////
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
			//else if (lined > distance * 1000)// 超出，认为不可通信
			//{
			//	Fvalue[0][j] = 1000;	// 超出通信范围，认为损耗极大，设为衰减1000dB
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
	xt=Setparameter[0]; // 发射天线经度
	yt=Setparameter[1]; // 发射天线纬度
	Gt=Setparameter[2]; // 发射天线最大天线增益
	ht=Setparameter[3]; // 发射天线高度

	xr=Receparameter[0]; // 接收天线经度
	yr=Receparameter[1]; // 接收天线纬度
	Gr=Receparameter[2]; // 接收天线最大天线增益
	hr=Receparameter[3]; // 接收天线高度

	row = dimension[0];    // 地理信息数据的行数
	column = dimension[1]; // 地理信息数据的列数
	acc = IPara.pSize;	   // 预测颗粒度
	distance = IPara.Elength;// 有效预测范围
	if(IPara.aMark == 0)// 淡水
	{
		aa=0.001;		// 有效电导率
		e=80;			// 相对介电常数
	}
	if(IPara.aMark == 1)// 海水
	{
		aa=4.0;			// 有效电导率
		e=80;			// 相对介电常数
	}
	if(IPara.aMark == 2)// 湿地
	{
		aa=0.008;		// 有效电导率
		e=10;			// 相对介电常数
	}
	if(IPara.aMark == 3)// 干地
	{
		aa=0.001;		// 有效电导率
		e=2;			// 相对介电常数
	}
	double **Fvalue; // 定义返回值二维指针
	double **GRID;   // 定义高程矩阵指针
	double dx,dy,dx1,dy1,xt1,yt1,xr1,yr1;
	int y3=1,y4=0;		// 极化参数和地面参数
	double re=8500000;  // 地球等效半径
	dy = fabs(Eledata[0][0].ycor-Eledata[1][1].ycor); // 相邻两个采样点的纬度差值
	dx = fabs(Eledata[0][0].xcor-Eledata[1][1].xcor); // 相邻两个采样点的经度差值
	// ======================== 地图中相邻采样点的水平间距和垂直间距m ============================ //
	//dx1 = dx * 6371000 * cos(Eledata[0][1].ycor/180*3.1415926) * 3.1415926/180.0; // 地图水平两点间距m
	//dy1 = dy * 6371000 * 3.1415926/180.0;                                         // 地图垂直两点间距m
    dx1 = 73.5624;			  // 地图水平两点间距m
	dy1 = 92.6043;			  // 地图垂直两点间距m
	// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = //
	xt1 = floor(fabs(xt - Eledata[0][0].xcor)/dx);    // 发射点在地图矩阵上的位置
	yt1 = floor(fabs(yt - Eledata[0][0].ycor)/dy);
	xr1 = floor(fabs(xr - Eledata[0][0].xcor)/dx);    // 接收点在地图矩阵上的位置
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
	

		// 动态申请内存
		Fvalue = new double*[row];
		for(i=0;i < row;i++)
		{
			Fvalue[i] = new double[column]; 
		}
		//////////////////////// ITU-rp526场强预测模型类 ////////////////////////////////////////////
		ITUrp RaosheITU; 
		double lined;
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < column; j++)
			{
				lined = sqrt(((xt1 - i) * dx1) * ((xt1 - i) * dx1) + ((yt1 - j) * dy1) * ((yt1 - j) * dy1)); //发射站到地图某点的距离
				if (lined == 0)
				{
					Fvalue[i][j] = 0.0;
				}
				else if (lined > distance*1000)// 超出，认为不可通信
				{
					Fvalue[i][j] = 1000;	// 超出通信范围，认为损耗极大，设为衰减1000dB
				}
				else if (lined < 1000)
				{
					Fvalue[i][j] = 32.45 + 20 * log10(lined / 1000) + 26 * log10(f); //自由空间传播
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
	xt=Setparameter[0]; // 发射天线经度
	yt=Setparameter[1]; // 发射天线纬度
	Gt=Setparameter[2]; // 发射天线最大天线增益
	ht=Setparameter[3]; // 发射天线高度
	Azimutht = Setparameter[4]/180*PI; // 发射天线方位角

	xr=Receparameter[0]; // 接收天线经度
	yr=Receparameter[1]; // 接收天线纬度
	Gr=Receparameter[2]; // 接收天线最大天线增益
	hr=Receparameter[3]; // 接收天线高度
	Azimuthr = Receparameter[4]/180*PI;// 接收天线方位角

	row = dimension[0];    // 地理信息数据的行数
	column = dimension[1]; // 地理信息数据的列数
	acc = Para.pSize*1000; // 预测颗粒度 单位m
	distance = Para.Elength;// 有效预测范围
	zetaN = Para.deterN;	// 大气折射率指数梯度
	PW = Para.pw;			// 最差月份的时间百分比
	Diameter = Para.diameter;// 天线直径
	double fre = f/1000;
	double *Fvalue; // 定义返回值一维指针
	double **GRID;   // 定义高程矩阵指针
	double dx,dy,dx1,dy1,xt1,yt1,xr1,yr1,Gaint;
	double re=8500000;  // 地球等效半径
	double fai = (yt+yr)/2.0;
	dy = fabs(Eledata[0][0].ycor-Eledata[1][1].ycor); // 相邻两个采样点的纬度差值
	dx = fabs(Eledata[0][0].xcor-Eledata[1][1].xcor); // 相邻两个采样点的经度差值
	// ======================== 地图中相邻采样点的水平间距和垂直间距m ============================ //
	//dx1 = dx * 6371000 * cos(Eledata[0][1].ycor/180*3.1415926) * 3.1415926/180.0; // 地图水平两点间距m
	//dy1 = dy * 6371000 * 3.1415926/180.0;                                         // 地图垂直两点间距m
    dx1 = 73.5624;			  // 地图水平两点间距m
	dy1 = 92.6043;			  // 地图垂直两点间距m
	// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = //
	xt1 = floor(fabs(xt - Eledata[0][0].xcor)/dx);    // 发射点在地图矩阵上的位置
	yt1 = floor(fabs(yt - Eledata[0][0].ycor)/dy);
	xr1 = floor(fabs(xr - Eledata[0][0].xcor)/dx);    // 接收点在地图矩阵上的位置
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
	// = = = = = = = = = = = = = = = = = = = = 微波链路衰减预测模型 = = = = = = = = = = = = = = = = = = //		
	double lamde = 3.0/fre;
	//Gaint = 10*log10(0.6*(PI*Diameter/lamde)*(PI*Diameter/lamde));
	Gaint = cauclateGain(xt,yt,xr,yr,Azimutht,Diameter,fre);
	double Length_L,lined,spacX,spacY,xi,yj;
	int N_Sample;
	Length_L = sqrt(((xt1 - xr1) * dx1) * ((xt1 - xr1) * dx1) + ((yt1 - yr1) * dy1) * ((yt1 - yr1) * dy1));// 通信链路长度
	N_Sample = (int)floor(Length_L/acc);		// 链路总的计算点数
	spacX = ((xt1 - xr1) * dx1)/N_Sample;	// 相邻两个计算点的经度差值
	spacY = ((yt1 - yr1) * dy1)/N_Sample;	// 相邻两个计算点的纬度差值
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
		else if (lined > distance*1000)// 超出，认为不可通信
		{
			Fvalue[j] = 1000;	// 超出通信范围，认为损耗极大，设为衰减1000dB
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
	xt=Setparameter[0]; // 发射天线经度
	yt=Setparameter[1]; // 发射天线纬度
	Gt=Setparameter[2]; // 发射天线最大天线增益
	ht=Setparameter[3]; // 发射天线高度

	xr=Receparameter[0]; // 接收天线经度
	yr=Receparameter[1]; // 接收天线纬度
	Gr=Receparameter[2]; // 接收天线最大天线增益
	hr=Receparameter[3]; // 接收天线高度

	row = dimension[0];    // 地理信息数据的行数
	column = dimension[1]; // 地理信息数据的列数
	acc = IPara.pSize;	   // 预测颗粒度
	distance = IPara.Elength;// 有效预测范围
	if(IPara.aMark == 0)// 淡水
	{
		aa=0.001;		// 有效电导率
		e=80;			// 相对介电常数
	}
	if(IPara.aMark == 1)// 海水
	{
		aa=4.0;			// 有效电导率
		e=80;			// 相对介电常数
	}
	if(IPara.aMark == 2)// 湿地
	{
		aa=0.008;		// 有效电导率
		e=10;			// 相对介电常数
	}
	if(IPara.aMark == 3)// 干地
	{
		aa=0.001;		// 有效电导率
		e=2;			// 相对介电常数
	}
	double ret_Value; // 定义返回值
	double **GRID;    // 定义高程矩阵指针
	double dx,dy,dx1,dy1,xt1,yt1,xr1,yr1;
	int y3=1,y4=0;		// 极化参数和地面参数
	double re=8500000;  // 地球等效半径
	dy = fabs(Eledata[0][0].ycor-Eledata[1][1].ycor); // 相邻两个采样点的纬度差值
	dx = fabs(Eledata[0][0].xcor-Eledata[1][1].xcor); // 相邻两个采样点的经度差值
	// ======================== 地图中相邻采样点的水平间距和垂直间距m ============================ //
	//dx1 = dx * 6371000 * cos(Eledata[0][1].ycor/180*3.1415926) * 3.1415926/180.0; // 地图水平两点间距m
	//dy1 = dy * 6371000 * 3.1415926/180.0;                                         // 地图垂直两点间距m
    dx1 = 73.5624;			  // 地图水平两点间距m
	dy1 = 92.6043;			  // 地图垂直两点间距m
	// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = //
	xt1 = floor(fabs(xt - Eledata[0][0].xcor)/dx);    // 发射点在地图矩阵上的位置
	yt1 = floor(fabs(yt - Eledata[0][0].ycor)/dy);
	xr1 = floor(fabs(xr - Eledata[0][0].xcor)/dx);    // 接收点在地图矩阵上的位置
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
	// = = = = = = = = = = = = = = = = = = = = ITU-rp526模型衰减预测 = = = = = = = = = = = = = = = = = = //
	double lined;
	lined = sqrt(((xt1 - xr1) * dx1) * ((xt1 - xr1) * dx1) + ((yt1 - yr1) * dy1) * ((yt1 - yr1) * dy1));// 通信链路长度
	ITUrp RaosheITU;// ITU-rp526场强预测模型类 
	if (lined == 0)
	{
		ret_Value = 0.0;
	}
	else if (lined > distance*1000)// 超出，认为不可通信
	{
		ret_Value = 1000;	// 超出通信范围，认为损耗极大，设为衰减1000dB
	}
	else if (lined < 1000)
	{
		ret_Value = 32.45 + 20 * log10(lined / 1000) + 26 * log10(f);
	}
	else
	{
		ret_Value = RaosheITU.ITUraoshe(xt1,yt1,xr1,yr1,GRID,row,column,acc,f,ht,hr,y3,y4,aa, e, re, dx1, dy1);
	}
		
	// = = = = = = = = = = = = = = = = = = = = 释放内存 = = = = = = = = = = = = = = = = = = = = = = = = = //	
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
	xt=Setparameter[0]; // 发射天线经度
	yt=Setparameter[1]; // 发射天线纬度
	Gt=Setparameter[2]; // 发射天线最大天线增益
	ht=Setparameter[3]; // 发射天线高度
	Azimutht = Setparameter[4]/180*PI; // 发射天线方位角

	xr=Receparameter[0]; // 接收天线经度
	yr=Receparameter[1]; // 接收天线纬度
	Gr=Receparameter[2]; // 接收天线最大天线增益
	hr=Receparameter[3]; // 接收天线高度
	Azimuthr = Receparameter[4]/180*PI;// 接收天线方位角

	row = dimension[0];    // 地理信息数据的行数
	column = dimension[1]; // 地理信息数据的列数
	acc = Para.pSize*1000; // 预测颗粒度 单位m
	distance = Para.Elength;// 有效预测范围
	zetaN = Para.deterN;	// 大气折射率指数梯度
	PW = Para.pw;			// 最差月份的时间百分比
	Diameter = Para.diameter;// 天线直径
	double fre = f/1000;
	double ret_Value; // 定义返回值一维指针
	double **GRID;   // 定义高程矩阵指针
	double dx,dy,dx1,dy1,xt1,yt1,xr1,yr1,Gaint,Gainr;
	double re=8500000;  // 地球等效半径
	double height1 = ht; 
	double fai = (yt+yr) / 2.0;
	dy = fabs(Eledata[0][0].ycor-Eledata[1][1].ycor); // 相邻两个采样点的纬度差值
	dx = fabs(Eledata[0][0].xcor-Eledata[1][1].xcor); // 相邻两个采样点的经度差值
	// ======================== 地图中相邻采样点的水平间距和垂直间距m ============================ //
	//dx1 = dx * 6371000 * cos(Eledata[0][1].ycor/180*3.1415926) * 3.1415926/180.0; // 地图水平两点间距m
	//dy1 = dy * 6371000 * 3.1415926/180.0;                                         // 地图垂直两点间距m
    dx1 = 73.5624;			  // 地图水平两点间距m
	dy1 = 92.6043;			  // 地图垂直两点间距m
	// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = //
	xt1 = floor(fabs(xt - Eledata[0][0].xcor)/dx);    // 发射点在地图矩阵上的位置
	yt1 = floor(fabs(yt - Eledata[0][0].ycor)/dy);
	xr1 = floor(fabs(xr - Eledata[0][0].xcor)/dx);    // 接收点在地图矩阵上的位置
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
	// = = = = = = = = = = = = = = = = = = = = 微波链路衰减预测模型 = = = = = = = = = = = = = = = = = = //		
	double lamde = 3.0/fre;
	//Gaint = 10*log10(0.6*(PI*Diameter/lamde)*(PI*Diameter/lamde));
	Gaint = cauclateGain(xt,yt,xr,yr,Azimutht,Diameter,fre);
	Gainr = cauclateGain(xr,yr,xt,yt,Azimuthr,Diameter,fre);
	double lined;
	lined = sqrt(((xt1 - xr1) * dx1) * ((xt1 - xr1) * dx1) + ((yt1 - yr1) * dy1) * ((yt1 - yr1) * dy1));// 通信链路长度
	
	attenuation_ra raITU;
	if (lined == 0)
	{
		ret_Value = 0.0;
	}
	else if (lined > distance*1000)// 超出，认为不可通信
	{
		ret_Value = 1000;	// 超出通信范围，认为损耗极大，设为衰减1000dB
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
