#ifndef DECAYFCDLL_H
#define DECAYFCDLL_H

struct Data
{
	double xcor;		// 经度
	double ycor;		// 纬度
	double eleva;		// 高程值
};
struct Iparameter
{
	double pSize;		// 预测颗粒度
	double Elength;		// 有效预测范围
	int aMark;			// 电导率和相对介电常数
};

struct FlagData
{
	double xcor;		// 经度
	double ycor;		// 纬度
	int flag;		    // 地形标志
};

struct ParWave
{
	double pSize;		// 预测颗粒度
	double Elength;		// 有效预测范围
	double deterN;		// 大气折射率指数梯度
	double pw;			// 最差月份的时间百分比
	double diameter;	// 天线直径
};

extern "C" double _declspec(dllexport) **DecayDiffraction(double[] ,double[],Data **,int[],double,Iparameter,double,double,int);

//extern "C" double _declspec(dllexport) *DecayMicrowave(double[],double[],Data **,int[],double,ParWave);
//
//extern "C" double _declspec(dllexport) PointDiffraction(double[],double[],Data **,int[],double,Iparameter);
//
//extern "C" double _declspec(dllexport) PointMicrowave(double[],double[],Data **,int[],double,ParWave);	
									   
#endif