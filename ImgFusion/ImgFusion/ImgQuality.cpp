#include "Function.h"
#include "Macro.h"
#include "ImgQuality.h"
#include <Windows.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <atlstr.h>  
#include <math.h>

using namespace std;
/***************************************************�ҶȻ�����Ϣ***********************************************************************************/
double CalAverage(unsigned char *ptImage,int width,int height) 
{
	if (ptImage == NULL)
	{
		MessageBox(NULL,_T("Ӱ��ָ��ָ��ʧ�ܣ����飡"),_T("ϵͳ��ʾ��"),MB_OK|MB_ICONERROR);
		return false;
	}
	LONG64 Total = 0;
	double ImgAver = 0.0;
	int l = 0;
	for (int i = 0;i < height;i++)
	{
		for (int j = 0;j < width;j++)
		{
			Total += ptImage[i * width + j];
			l++;
		}
	}
	//cout << (int) ptImage[100 * PanWidth + 100] << endl;
	//cout<<Total<< "\t" <<l<<endl;
	ImgAver = static_cast<double>(Total)/(width * height);//��ƽ��
	//�ؼ���static_cast ������������ͬ����֮�����ǿ��ת��,����û������ʱ����.
	printf("Ӱ���ֵΪ��%.6g\n",ImgAver);
	return ImgAver;
}

double CalStandardDevi(unsigned char *ptImage,int width,int height)
{
	if (ptImage == NULL)
	{
		MessageBox(NULL,_T("Ӱ��ָ��ָ��ʧ�ܣ����飡"),_T("ϵͳ��ʾ��"),MB_OK|MB_ICONERROR);
		return false;
	}
	double aver = CalAverage(ptImage,width,height);
	double Var = 0;
	double SD = 0;
	for (int i = 0;i < height;i++)
	{
		for (int j = 0;j < width;j++)
		{
			Var += pow((ptImage[i * width + j] - aver),2);
		}
	}
	int SumCount = width*height;
	SD = sqrt(Var/(double)SumCount);
	printf("Ӱ���׼��Ϊ��%.5g\n",SD);
	return SD;
}
/***************************************************������***********************************************************************************/
double CalSpatialFreq(unsigned char *ptImage,int width,int height)
{
	if (ptImage == NULL)
	{
		MessageBox(NULL,_T("Ӱ��ָ��ָ��ʧ�ܣ����飡"),_T("ϵͳ��ʾ��"),MB_OK|MB_ICONERROR);
		return false;
	}
	double SumRow = 0.0;
	double SumColumn = 0.0;
	double SF = 0.0;
	for (int i = 0;i < height;i++)
	{
		for (int j = 1;j < width;j++)
		{
			SumRow += pow((double)(ptImage[i * width + j] - ptImage[i*width + j -1]),2);
		}
	}
	for (int m = 1;m < height;m++)
	{
		for (int n = 0;n < width;n++)
		{
			SumColumn += pow((double)(ptImage[m * width + n] - ptImage[(m - 1)*width + n]),2);
		}
	}
	int SumCount = width * height;
	double RF = sqrt(SumRow/(double)SumCount);
	double CF = sqrt(SumColumn/(double)SumCount);
    SF = sqrt(pow(SF,2) + pow(CF,2));
	printf("Ӱ��ռ�Ƶ��Ϊ��%.4g\n",SF);
	return SF;
}

double CalAverGrad(unsigned char *ptImage,int width,int height)
{
	if (ptImage == NULL)
	{
		MessageBox(NULL,_T("Ӱ��ָ��ָ��ʧ�ܣ����飡"),_T("ϵͳ��ʾ��"),MB_OK|MB_ICONERROR);
		return false;
	}
	double AG = 0;
	double dx = 0;
	double dy = 0;
	double ds = 0;
	for (int i = 0;i < height - 1; i++)
	{
		for (int j = 0;j < width - 1; j++)
		{
			dx = ptImage[i * width + j + 1] - ptImage[i * width + j];
			dy = ptImage[(i + 1) * width + j] - ptImage[i * width + j];
			ds += sqrt((pow(dx,2) + pow(dy,2))/2);
		}
	}
	AG = ds/((width - 1)*(height - 1));


	//double AG = 0.0;
	//double AGSum = 0.0;
	//double AGx = 0.0;
	//double AGy = 0.0;
	//
	//unsigned char *temp = new unsigned char[FilterWidth * FilterWidth];
	//memset(temp,0,sizeof(unsigned char)*FilterWidth*FilterWidth);
 //   //3*3ģ�� Soble����
	//int xFilter[FilterWidth * FilterWidth] = { -1,-2,-1,
	//	                                        0, 0, 0,
	//                                          	1, 2, 1
	//                                         };
	//int yFilter[FilterWidth * FilterWidth] = { -1,0,1,
	//	                                       -2,0,2,
	//	                                       -1,0,1
	//                                         };
	//
	//for (int i = FilterWidth/2;i < PanHeight - FilterWidth/2;++i)       
	//{
	//	for (int j = FilterWidth/2;j < PanWidth - FilterWidth/2;++j)
	//	{
	//		int index=0;
	//		for (int k = -FilterWidth/2;k <= FilterWidth/2;++k)
	//		{
	//			for (int h = -FilterWidth/2;h <= FilterWidth/2;++h)
	//			{
	//				temp[index] = ptImage[(i + k) * PanWidth + (j + h)] ;
	//				index++;
	//			}
	//		}
	//		for (int m = 0;m < FilterWidth;m++)
	//		{
	//			for (int n = 0;n < FilterWidth;n++)
	//			{
	//				AGx += temp[m * FilterWidth + n] * xFilter[m * FilterWidth + n];
	//				AGy += temp[m * FilterWidth + n] * yFilter[m * FilterWidth + n];
	//				
	//			}
	//		}
	//		AGSum += sqrt((pow(AGx,2) + pow(AGy,2))/2);
	//	}
	//}
	//AG = AGSum/(double)((PanWidth - 1)*(PanHeight - 1));
	//delete []temp;
	printf("Ӱ��ƽ���ݶ�Ϊ��%.4g\n",AG);
	
	return AG;
}
/***************************************************��Ϣ��***********************************************************************************/
double CalEntropy(unsigned char *ptImage,int width,int height)
{
	if (ptImage == NULL)
	{
		MessageBox(NULL,_T("Ӱ��ָ��ָ��ʧ�ܣ����飡"),_T("ϵͳ��ʾ��"),MB_OK|MB_ICONERROR);
		return false;
	}
	long lCount[ColorLen];//longָ��long int����
	
	double SumTemp = 0.0;
	double Entropy = 0.0;
	memset(lCount,0,sizeof(long)*ColorLen);
	
	for (int i = 0;i < height;i++)
	{
		for (int j = 0;j < width;j++)
		{
			lCount[ptImage[i * width + j]]++;
		}
	}
	int SumCount = width * height;
	double percent = 0.0;
	for (int l = 0;l < ColorLen;l++)
	{
		percent = lCount[l] /(double)SumCount;                          //����������ת������Ȼ�ض�Ӱ����
		//cout << percent << endl;
		SumTemp += percent*(logl(percent + 0.000000001)/logl(2));     //����һ�����׹�ʽ,Ϊ��ֹ����Ϊ0��һ����Сֵ
	}
	Entropy = -SumTemp;
	printf("Ӱ����Ϣ��Ϊ��%.4g\n",Entropy);
	return Entropy;
}
/***************************************************���ױ����***********************************************************************************/
double CalCorrelationCoef(unsigned char *dstImage,unsigned char *refImage,int width,int height)
{
	double uf = CalAverage(dstImage,width,height);
	double ur = CalAverage(refImage,width,height);
	double tempA = 0.0;
	double tempB = 0.0;
	double tempC = 0.0;
	double tempD = 0.0;
	double CC = 0.0;
	for (int i = 0;i < height;i++)
	{
		for (int j = 0;j < width;j++)
		{
			tempA += (dstImage[i * width + j] - uf) * (refImage[i * width + j] - ur);
			tempB += pow((dstImage[i * width + j] - uf),2);
			tempC += pow((refImage[i * width + j] - ur),2);
		}
	}
	tempD = sqrt(tempB * tempC);
	CC = tempA/tempD;
	printf("Ӱ�������ϵ��Ϊ��%.5g \n",CC);
	return CC;
}


double CalDeviationIndex(unsigned char *dstImage,unsigned char *refImage,int width,int height)
{
	double temp = 0.0;
	double D = 0.0;
	for (int i = 0;i < height;i++)
	{
		for (int j = 0;j < width;j++)
		{
			temp += abs((dstImage[i * width + j] /refImage[i * width + j])-1);
		}
	}
	D = temp/((width - 1)*(height - 1));
	printf("Ӱ����ƫ��ָ��Ϊ��%.5g \n",D);
	return D;
}

double CalDistortion(unsigned char *dstImage,unsigned char *refImage,int width,int height)
{

	double temp = 0.0;
	double DD = 0.0;
	for (int i = 0;i < height;i++)
	{
		for (int j = 0;j < width;j++)
		{
			temp += abs(dstImage[i * width + j] - refImage[i * width + j]);
		}
	}
	DD = temp/((width - 1)*(height - 1));
	printf("Ӱ����Ť����Ϊ��%.5g \n",DD);
	return DD;
}


















