#ifndef DECAYFCDLL_H
#define DECAYFCDLL_H

struct Data
{
	double xcor;		// ����
	double ycor;		// γ��
	double eleva;		// �߳�ֵ
};
struct Iparameter
{
	double pSize;		// Ԥ�������
	double Elength;		// ��ЧԤ�ⷶΧ
	int aMark;			// �絼�ʺ���Խ�糣��
};

struct FlagData
{
	double xcor;		// ����
	double ycor;		// γ��
	int flag;		    // ���α�־
};

struct ParWave
{
	double pSize;		// Ԥ�������
	double Elength;		// ��ЧԤ�ⷶΧ
	double deterN;		// ����������ָ���ݶ�
	double pw;			// ����·ݵ�ʱ��ٷֱ�
	double diameter;	// ����ֱ��
};

extern "C" double _declspec(dllexport) **DecayDiffraction(double[] ,double[],Data **,int[],double,Iparameter,double,double,int);

//extern "C" double _declspec(dllexport) *DecayMicrowave(double[],double[],Data **,int[],double,ParWave);
//
//extern "C" double _declspec(dllexport) PointDiffraction(double[],double[],Data **,int[],double,Iparameter);
//
//extern "C" double _declspec(dllexport) PointMicrowave(double[],double[],Data **,int[],double,ParWave);	
									   
#endif