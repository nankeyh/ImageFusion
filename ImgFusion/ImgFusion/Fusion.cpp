#include "Fusion.h"
#include <vector>
#include <atlstr.h>
#include "windows.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <time.h>
#include <string>
#include <assert.h>
#include <math.h>
#include <list>

using namespace std;

//////////////////////////////////////////////////////////////////////////
/*                              HPF（高通滤波）                         */
//////////////////////////////////////////////////////////////////////////
void HPFFusion(const char* msfilepath,const char* panfilepath,char* savefusionpath)  
{  
	int width =0;
	int height = 0;
	int nRastercount=0;
	//读入多光谱影像
	unsigned char **pMSImage = GDALreadImage(width,height,nRastercount,msfilepath);
	//滤波器的选择3×3
	int FilterWidth = 3;
	//读入全色影像
	unsigned char *pPANImage = OpenTif(panfilepath);
	//中间过程
	double *temppan = new double[width*height];
	memset(temppan,0,sizeof(double)*width*height);
	//滤波器
	unsigned char *temp = new unsigned char[FilterWidth*FilterWidth];
	memset(temp,0,sizeof(unsigned char)*FilterWidth*FilterWidth);

	//i,j从1开始到width-1结束，这样i*width+j才能表示滤波器中间值。相当于外围一圈继续用原值，内部经过滤波重新赋值
	for (int i = FilterWidth/2;i < height - FilterWidth/2;++i)       
	{
		for (int j = FilterWidth/2;j < width - FilterWidth/2;++j)
		{
			int index=0;double sum = 0;
			for (int k = -FilterWidth/2;k <= FilterWidth/2;++k)
			{
				for (int h = -FilterWidth/2;h <= FilterWidth/2;++h)
				{
					temp[index] = pPANImage[(i + k) * width + (j + h)];//遍历获得滤波器中各个灰度值
					index++;
					sum += temp[index];
				}
			}
			temppan[i * width + j] = sum - FilterWidth*FilterWidth*temp[FilterWidth*FilterWidth/2];
		}
	}
	delete []temp; temp = NULL;
	delete []pPANImage;pPANImage = NULL;


	double **tempfusion = new double *[nRastercount];
	unsigned char **FusionImg = new unsigned char *[nRastercount];
	for (int k = 0;k < nRastercount;k++)
	{
		tempfusion[k] = new double[width*height];
		memset(tempfusion[k],0,sizeof(double)*width*height);
		FusionImg[k] = new unsigned char[width*height];
		memset(FusionImg[k],0,sizeof(unsigned char)*width*height);
	}
	//把高频成分加到多光谱图像上面
	for(int i = 0;i < height;i++)
	{
		for(int j = 0;j < width;j++)
		{	
			for (int n = 0;n < nRastercount;n++)
			{
				//tempfusion[n][i*width+j] = (double)pMSImage[n][i*width+j] ;
				tempfusion[n][i*width+j] = (double)pMSImage[n][i*width+j] + temppan[i*width+j];
			}
			//tempfusion[0][i*width+j] += temppan[i*width+j];
			//只在某一个波段上加高频信息，结果影像是模糊的（波段间分辨率不一致），不能达到增强影像分辨率的要求
		}
	}
	delete []temppan;temppan = NULL;
	for (int i = 0;i < nRastercount;i++)
	{
		delete []pMSImage[i];
		pMSImage[i] = NULL;
	}
	delete []pMSImage;pMSImage = NULL;
	linearStretch(tempfusion,width,height,nRastercount);
	//同时加高频信息再线性拉伸，各波段的像素值基本持平，从而使得融合影像白化并致使光谱信息丢失严重
	for (int k = 0;k < nRastercount;k++)
	{
		FusionImg[k] = new unsigned char[width*height];
		for (int s = 0;s < width*height;s++)
		{
			FusionImg[k][s] = (unsigned char)tempfusion[k][s];
		}
	}

	for (int k = 0; k < nRastercount;k++)
	{
		delete []tempfusion[k];
		tempfusion[k] = NULL;
	}
	delete []tempfusion;
	tempfusion = NULL;

	SaveImage(FusionImg,savefusionpath,width,height,nRastercount);

	for (int k = 0; k < nRastercount;k++)
	{
		delete []FusionImg[k];
		FusionImg[k] = NULL;
	}
	delete []FusionImg;
	FusionImg = NULL;
	return ;  	
}

//////////////////////////////////////////////////////////////////////////
/*                    乘积性变换融合算法（Multiply以及改进算法）        */
//////////////////////////////////////////////////////////////////////////
void MultiplyFusion(const char* msfilepath,const char* panfilepath,char* savefusionpath,int flag)
{
	int width =0;
	int height = 0;
	int nRastercount=0;
	double RGBsum = 0.0;
	double panvalue = 0.0;
	//读入多光谱影像
	unsigned char **pMSImage = GDALreadImage(width,height,nRastercount,msfilepath);
	//读入全色影像
	unsigned char* pPANImage = OpenTif(panfilepath);
	//中间过程
	double**tempimg = new double*[nRastercount];
	for (int k = 0;k < nRastercount;k++)
	{
		tempimg[k] = new double[width*height];
		memset(tempimg[k],0,sizeof(double)*width*height);
	}
	switch(flag)
	{
	case 0:
		cout<<"乘积性变换"<<endl;	
		for (int i = 0;i < height;i++)
		{
			for (int j = 0;j < width;j++)
			{
				for (int k = 0;k < nRastercount;k++)
				{
					tempimg[k][i*width+j] = pPANImage[i*width+j]*pMSImage[k][i*width+j];
				}
			}
		}
		linearStretch(tempimg,width,height,nRastercount);
		break;
	case 1:
		cout<<"线性复合与乘积变换"<<endl;
		for (int i = 0;i < height;i++)
		{
			for (int j = 0;j < width;j++)
			{
				panvalue += pPANImage[i*width+j];
				for (int k = 0;k < nRastercount;k++)
				{
					RGBsum += pMSImage[k][i*width+j];
				}
				
			}
		}
		if (panvalue == 0.0)panvalue = 1.0;
		if (RGBsum == 0.0)RGBsum = 1.0;
		for (int i = 0;i < height;i++)
		{
			for (int j = 0;j < width;j++)
			{
				double m = (RGBsum-3*panvalue)/(width*height);
				for (int k = 0;k < nRastercount;k++)
				{
					tempimg[k][i*width+j] = m+(sqrt((double)(pMSImage[k][i*width+j]*pPANImage[i*width+j])))*0.95;
				}
			}
		}
		linearStretch(tempimg,width,height,nRastercount);
		break;
	default:
		cout<<"flag does not exist."<<endl;
		break;
	}

	for (int k = 0; k < nRastercount;k++)
	{
		delete []pMSImage[k];
		pMSImage[k] = NULL;
	}
	delete []pMSImage;
	pMSImage = NULL;

	delete []pPANImage;
	pPANImage = NULL;
	//linearStretch(tempimg,width,height,nRastercount);

	unsigned char **FusionImg = new unsigned char *[nRastercount];
	for (int k = 0;k < nRastercount;k++)
	{
		FusionImg[k] = new unsigned char[width*height];
		for (int s = 0;s < width*height;s++)
		{
			FusionImg[k][s] = (unsigned char)tempimg[k][s];
		}
	}

	for (int k = 0; k < nRastercount;k++)
	{
		delete []tempimg[k];
		tempimg[k] = NULL;
	}
	delete []tempimg;
	tempimg = NULL;
	SaveImage(FusionImg,savefusionpath,width,height,nRastercount);

	for (int k = 0; k < nRastercount;k++)
	{
		delete []FusionImg[k];
		FusionImg[k] = NULL;
	}
	delete []FusionImg;
	FusionImg = NULL;

	return ;
}

//////////////////////////////////////////////////////////////////////////
/*                   彩色变换融合算法（IHS融合方法以及相关算法）        */
//////////////////////////////////////////////////////////////////////////
void IHSFusion(const char* msfilepath,const char* panfilepath,char* savefusionpath,int flag)
{
	int width =0;
	int height = 0;
	int nRastercount=0;
	//读入多光谱影像
	unsigned char **pMSImage = GDALreadImage(width,height,nRastercount,msfilepath);
	if (nRastercount != 3)
	{
		cout<<"读入的低分辨率影像波段数有误，不能进行该融合操作！"<<endl;
		for (int k = 0; k < nRastercount;k++)
		{
			delete []pMSImage[k];
			pMSImage[k] = NULL;
		}
		delete []pMSImage;
		pMSImage = NULL;
		return;
	}
	//读入全色影像
	unsigned char* pPANImage = OpenTif(panfilepath);
	//中间过程
	double**tempimg = new double*[nRastercount];//加入中间过程，结果影像的效果变化不大？？？
	for (int k = 0;k < nRastercount;k++)
	{
		tempimg[k] = new double[width*height];
		memset(tempimg[k],0,sizeof(double)*width*height);
	}
	double RGBsum = 0.0;
	double dmax = 0.0;
	double dmin = 255;
	double C = 0,PI = 3.1415926;
	double v1,v2,fI,fH,fS,br,bg,bb;
	double*dRGB = new double[nRastercount];
	double*dRGBt = new double[nRastercount];
	switch(flag)
	{
	case 0:
		cout<<"IHS变换"<<endl;	
		for (int i = 0;i < height;i++)
		{
			for (int j = 0;j < width;j++)
			{
				memset(dRGB,0,sizeof(double)*nRastercount);
				memset(dRGBt,0,sizeof(double)*nRastercount);
				for (int k = 0;k < nRastercount;k++)
				{
					dRGB[k] = pMSImage[k][i*width+j];
					RGBsum += dRGB[k];
				}
				if (RGBsum == 0.0)
				{
					RGBsum = 1.0;
				}
				//进行正变换
				fI = RGBsum/3;
				max_min_value(dRGB,nRastercount,dmax,dmin);
				fS = (RGBsum - 3*dmin)/RGBsum;
				double angle = acos(0.5*(2*dRGB[0]-dRGB[1]-dRGB[2])/(sqrt((dRGB[0]-dRGB[1])*(dRGB[0]-dRGB[1])+(dRGB[0]-dRGB[2])*(dRGB[1]-dRGB[2]))));
				angle = angle*180.0/PI;
				if (dRGB[2]>dRGB[1])//BG比较
				{
					fH = (360-angle)/360;
				}
				else
				{
					fH = angle/360;
				}
				
				dRGBt[0] = fI;dRGBt[1] = fS;dRGBt[2] = fH;
				//用高空间分辨率的图像替换Intensity亮度分量
				dRGBt[0] = (double)pPANImage[i*width+j];
				//进行逆变换
				dRGBt[2] *= 360;
				if (dRGBt[2] < 0)
				{
					dRGBt[2] += 360;
				}
				if (dRGBt[2]>=0&&dRGBt[2]<120)
				{
					br = dRGBt[0]*(1+((dRGBt[1]*cos(dRGB[2]*PI/180))/cos((60-dRGB[2])*PI/180)));			
					bb = dRGBt[0]*(1-dRGBt[1]);
					bg = dRGBt[0]*3.0-(br+bb);
				}
				if (dRGBt[2]>=120&&dRGBt[2]<240)
				{
					dRGBt[2]-=120;
					br = dRGBt[0]*(1-dRGBt[1]);
					bg = dRGBt[0]*(1+((dRGBt[1]*cos(dRGB[2]*PI/180))/cos((60-dRGB[2])*PI/180)));			
					bb = dRGBt[0]*3.0-(br+bg);
				}
				if (dRGBt[2]>=240&&dRGBt[2]<=360)
				{
					dRGBt[2]-=240;
					bg = dRGBt[0]*(1-dRGBt[1]);
					bb = dRGBt[0]*(1+((dRGBt[1]*cos(dRGB[2]*PI/180))/cos((60-dRGB[2])*PI/180)));			
					br = dRGBt[0]*3.0-(br+bg);
				}
				if(br>255) br = 255;
				if(br<0)   br = 0;
				if(bg>255) bg = 255;
				if(bg<0)  bg = 0;
				if(bb>255) bb = 255;
				if(bb<0)  bb = 0;
				dRGB[0] = br;dRGB[1] = bg;dRGB[2] = bb;
				//放在内存中
				for (int k = 0;k < nRastercount;k++)
				{
					tempimg[k][i*width+j] = dRGB[k];
				}
				dmax = 0;dmin = 255;
				RGBsum = 0.0;
			}
		}
		break;
	case 1:
		cout<<"进行三角IHS变换"<<endl;	
		for (int i = 0;i < height;i++)
		{
			for (int j = 0;j < width;j++)
			{
				memset(dRGB,0,sizeof(double)*nRastercount);
				memset(dRGBt,0,sizeof(double)*nRastercount);
				for (int k = 0;k < nRastercount;k++)
				{
					dRGB[k] = pMSImage[k][i*width+j];
					RGBsum += dRGB[k];
				}
				if (RGBsum == 0.0)
				{
					RGBsum = 1.0;
				}
				//进行正变换
				fI = RGBsum/3;
				max_min_value(dRGB,nRastercount,dmax,dmin);
				fS = (RGBsum - 3*dmin)/RGBsum;
				if (dmin == dRGB[0])//R最小
				{
					//fS = (RGBsum - 3*dRGB[0])/RGBsum;
					if((RGBsum-3*dRGB[0]) != 0) 
					{
						fH = (dRGB[2]-dRGB[0])/(RGBsum-3*dRGB[0])+1;
					}
					else //加判断预防三数相等情况，一般都是true
					{
						fH = (dRGB[2]-dRGB[0])/(RGBsum-3*dRGB[0]+1)+1;	
					}
				}
				if (dmin == dRGB[2])//B最小
				{
					//fS = (RGBsum - 3*dRGB[2])/RGBsum;
					if((RGBsum-3*dRGB[2]) != 0) 
					{
						fH = (dRGB[1]-dRGB[2])/(RGBsum-3*dRGB[2]);
					}
					else 
					{
						fH = (dRGB[1]-dRGB[2])/(RGBsum-3*dRGB[2]+1);	
					}
				}
				if (dmin == dRGB[1])//G最小
				{
					//fS = (RGBsum - 3*dRGB[1])/RGBsum;
					if((RGBsum-3*dRGB[1]) != 0) 
					{
						fH = (dRGB[0]-dRGB[1])/(RGBsum-3*dRGB[1])+2;
					}
					else 
					{
						fH = (dRGB[0]-dRGB[1])+2;	
					}
				}
				dRGBt[0] = fI;dRGBt[1] = fS;dRGBt[2] = fH;
				//用高空间分辨率的图像替换Intensity亮度分量
				dRGBt[0] = (double)pPANImage[i*width+j];
				//进行逆变换
				if (dmin == dRGB[0])
				{
					br = dRGBt[0]*(1-dRGBt[1]);
					bg = dRGBt[0]*(1+5*dRGBt[1]-3*dRGBt[1]*dRGBt[2]);
					bb = dRGBt[0]*(1-4*dRGBt[1]+3*dRGBt[1]*dRGBt[2]);
				}				
				if (dmin == dRGB[2])
				{
					br = dRGBt[0]*(1+2*dRGBt[1]-3*dRGBt[1]*dRGBt[2]);
					bg = dRGBt[0]*(1-dRGBt[1]+3*dRGBt[1]*dRGBt[2]);
					bb = dRGBt[0]*(1-dRGBt[1]);
				}
				if (dmin == dRGB[1])
				{
					br = dRGBt[0]*(1-7*dRGBt[1]+3*dRGBt[1]*dRGBt[2]);
					bg = dRGBt[0]*(1-dRGBt[1]);
					bb = dRGBt[0]*(1+8*dRGBt[1]-3*dRGBt[1]*dRGBt[2]);
				}
				if(br>255) br = 255;
				if(br<0)   br = 0;
				if(bg>255) bg = 255;
				if(bg<0)  bg = 0;
				if(bb>255) bb = 255;
				if(bb<0)  bb = 0;
				dRGB[0] = br;dRGB[1] = bg;dRGB[2] = bb;
				//放在内存中
				for (int k = 0;k < nRastercount;k++)
				{
					tempimg[k][i*width+j] = dRGB[k];
				}
				dmax = 0;dmin = 255;
				RGBsum = 0.0;
			}
		}
		break;
	case 2:
		cout<<"进行普通圆柱体IHS变换(间接IHS)"<<endl;	
		for (int i = 0;i < height;i++)
		{
			for (int j = 0;j < width;j++)
			{
				memset(dRGB,0,sizeof(double)*nRastercount);
				memset(dRGBt,0,sizeof(double)*nRastercount);
				for (int k = 0;k < nRastercount;k++)
				{
					dRGB[k] = pMSImage[k][i*width+j];
				}
				//进行正变换//G、B比较
				if(dRGB[1]>=dRGB[2]) C = 0;
				if(dRGB[1]<dRGB[2]) C = PI;
				fI = (1/sqrt(3.0))*(dRGB[2]+dRGB[1]+dRGB[0]);
				v1 = (dRGB[2]+dRGB[1]-2*dRGB[0])/(sqrt(6.0));
				v2 = (dRGB[2]-dRGB[1])/(sqrt(2.0));
				fH = atan(v1/v2) + C;
				fS = sqrt(v1*v1+v2*v2);			   
				dRGBt[0] = fI;dRGBt[1] = fS;dRGBt[2] = fH;
				//用高空间分辨率的图像替换Intensity亮度分量
				dRGBt[0] = pPANImage[i*width+j];
				//进行逆变换
				br = (1/sqrt(3.0))*dRGBt[0]-(2/sqrt(6.0))*v1;
				bg = (1/sqrt(3.0))*dRGBt[0]+(1/sqrt(6.0))*v1-(1/sqrt(2.0))*v2;
				bb = (1/sqrt(3.0))*dRGBt[0]+(1/sqrt(6.0))*v1+(1/sqrt(2.0))*v2;		
				br *= 1.5;
				bg *= 1.5;
				bb *= 1.5;
				if(br>255) br = 255;
				if(br<0)   br = 0;
				if(bg>255) bg = 255;
				if(bg<0)  bg = 0;
				if(bb>255) bb = 255;
				if(bb<0)  bb = 0;
				dRGB[0] = br;dRGB[1] = bg;dRGB[2] = bb;
				//放在内存中
				for (int k = 0;k < nRastercount;k++)
				{
					tempimg[k][i*width+j] = dRGB[k];
				}
				dmax = 0;dmin = 255;
				RGBsum = 0.0;
			}
		}
		break;
	case 3:
		cout<<"进行拉伸的圆柱体IHS变换"<<endl;
		for (int i = 0;i < height;i++)
		{
			for (int j = 0;j < width;j++)
			{
				memset(dRGB,0,sizeof(double)*nRastercount);
				memset(dRGBt,0,sizeof(double)*nRastercount);
				for (int k = 0;k < nRastercount;k++)
				{
					dRGB[k] = pMSImage[k][i*width+j];
				}
				//进行正变换//G、B比较
				if(dRGB[1]>=dRGB[2]) C = 0;
				if(dRGB[1]<dRGB[2]) C = PI;//π
				fI = (1/sqrt(3.0))*(dRGB[2]+dRGB[1]+dRGB[0]);
				v1 = (dRGB[2]+dRGB[1]-2*dRGB[0])/(sqrt(6.0));
				v2 = (dRGB[2]-dRGB[1])/(sqrt(2.0));
				fH = atan(v1/v2) + C;
				fS = sqrt(v1*v1+v2*v2);			   
				dRGBt[0] = fI;dRGBt[1] = fS;dRGBt[2] = fH;
				//用高空间分辨率的图像替换Intensity亮度分量
				dRGBt[0] = pPANImage[i*width+j];
				//进行逆变换
				br = (1/sqrt(3.0))*dRGBt[0]-(2/sqrt(6.0))*v1;
				bg = (1/sqrt(3.0))*dRGBt[0]+(1/sqrt(6.0))*v1-(1/sqrt(2.0))*v2;
				bb = (1/sqrt(3.0))*dRGBt[0]+(1/sqrt(6.0))*v1+(1/sqrt(2.0))*v2;	
				/*br *= 1.5;
				bg *= 1.5;
				bb *= 1.5;*///????这样就拉伸了
				if(br>255) br = 255;
				if(br<0)   br = 0;
				if(bg>255) bg = 255;
				if(bg<0)  bg = 0;
				if(bb>255) bb = 255;
				if(bb<0)  bb = 0;
				dRGB[0] = br;dRGB[1] = bg;dRGB[2] = bb;
				//放在内存中
				for (int k = 0;k < nRastercount;k++)
				{
					tempimg[k][i*width+j] = dRGB[k];
				}
				dmax = 0;dmin = 255;
				RGBsum = 0.0;
			}
		}
		break;
	}
	delete []dRGB;dRGB=NULL;
	delete []dRGBt;dRGBt=NULL;

	for (int k = 0; k < nRastercount;k++)
	{
		delete []pMSImage[k];
		pMSImage[k] = NULL;
	}
	delete []pMSImage;
	pMSImage = NULL;

	delete []pPANImage;
	pPANImage = NULL;
	//linearStretch(tempimg,width,height,nRastercount);

	unsigned char **FusionImg = new unsigned char *[nRastercount];
	for (int k = 0;k < nRastercount;k++)
	{
		FusionImg[k] = new unsigned char[width*height];
		for (int s = 0;s < width*height;s++)
		{
			FusionImg[k][s] = (unsigned char)tempimg[k][s];
		}
	}

	for (int k = 0; k < nRastercount;k++)
	{
		delete []tempimg[k];
		tempimg[k] = NULL;
	}
	delete []tempimg;
	tempimg = NULL;
	SaveImage(FusionImg,savefusionpath,width,height,nRastercount);

	for (int k = 0; k < nRastercount;k++)
	{
		delete []FusionImg[k];
		FusionImg[k] = NULL;
	}
	delete []FusionImg;
	FusionImg = NULL;

	return ;
}

//求多个数的最大值和最小值
void max_min_value(double* array,int n, double& a,double& b)
{
	double dMax,dMin;
	dMax = dMin =array[0];
	double *p = array;
	for(int i=0; i<n;i++)
	{
		if(dMax < *p)
		{
			dMax = *p;
		}
		if(dMin > *p)
		{
			dMin = *p;
		}
		p++;
	}
	a=dMax;
	b=dMin;
}

//////////////////////////////////////////////////////////////////////////
/*                比值融合算法（Brovey融合方法以及相关改进算法）        */
//////////////////////////////////////////////////////////////////////////
void BroveyFusion(const char* msfilepath,const char* panfilepath,char* savefusionpath)
{
	int width =0;
	int height = 0;
	int nRastercount=0;
	//读入多光谱影像
	unsigned char **pMSImage = GDALreadImage(width,height,nRastercount,msfilepath);
	//读入全色影像
	unsigned char* pPANImage = OpenTif(panfilepath);
	//中间过程
	double**tempimg = new double*[nRastercount];//加入中间过程，结果影像的效果变化不大？？？
	for (int k = 0;k < nRastercount;k++)
	{
		tempimg[k] = new double[width*height];
		memset(tempimg[k],0,sizeof(double)*width*height);
	}
	double RGBsum = 0.0;
	for (int i = 0;i < height;i++)
	{
		for (int j = 0;j < width;j++)
		{
			for (int k = 0;k < nRastercount;k++)
			{
				RGBsum += pMSImage[k][i*width+j];
			}
			if (RGBsum == 0.0)
			{
				RGBsum = 1.0;
			}
			if (nRastercount == 3)
			{
				for (int k = 0;k < nRastercount;k++)
				{
					tempimg[k][i*width+j] = (double)((pPANImage[i*width+j]*pMSImage[k][i*width+j])/RGBsum);
				}
			} 
			else
			{
				for (int k = 0;k < nRastercount;k++)
				{
					tempimg[k][i*width+j] = (double)(sqrt((double)nRastercount)*(pPANImage[i*width+j]*pMSImage[k][i*width+j])/RGBsum);
				}
			}
			RGBsum = 0.0;
		}
	}
	for (int k = 0; k < nRastercount;k++)
	{
		delete []pMSImage[k];
		pMSImage[k] = NULL;
	}
	delete []pMSImage;
	pMSImage = NULL;

	delete []pPANImage;
	pPANImage = NULL;
	linearStretch(tempimg,width,height,nRastercount);

	unsigned char **FusionImg = new unsigned char *[nRastercount];
	for (int k = 0;k < nRastercount;k++)
	{
		FusionImg[k] = new unsigned char[width*height];
		for (int s = 0;s < width*height;s++)
		{
			FusionImg[k][s] = (unsigned char)tempimg[k][s];
		}
	}

	for (int k = 0; k < nRastercount;k++)
	{
		delete []tempimg[k];
		tempimg[k] = NULL;
	}
	delete []tempimg;
	tempimg = NULL;
	SaveImage(FusionImg,savefusionpath,width,height,nRastercount);

	for (int k = 0; k < nRastercount;k++)
	{
		delete []FusionImg[k];
		FusionImg[k] = NULL;
	}
	delete []FusionImg;
	FusionImg = NULL;

	return ;
}

//////////////////////////////////////////////////////////////////////////
/*                       基于相关系数的加权融合方法                     */
//////////////////////////////////////////////////////////////////////////
void CorrelationFusion(const char* msfilepath,const char* panfilepath,char* savefusionpath)
{
	int width =0;
	int height = 0;
	int nRastercount=0;
	//读入多光谱影像
	unsigned char **pMSImage = GDALreadImage(width,height,nRastercount,msfilepath);
	//读入全色影像
	unsigned char* pPANImage = OpenTif(panfilepath);
	//1.计算各波段均值
	double PanMean = 0;
	
	double *MSMean = new double[nRastercount];
	memset(MSMean,0,sizeof(double)*nRastercount);
	double *sumvalue = new double[nRastercount+1];
	memset(sumvalue,0,sizeof(double)*(nRastercount+1));
	for (int i = 0;i < height;i++)
	{
		for (int j = 0;j < width;j++)
		{
			sumvalue[nRastercount] += (double)pPANImage[i*width+j]; 
			for (int k = 0; k < nRastercount;k++)
			{
				sumvalue[k] = (double)pMSImage[k][i*width+j];
			}
		}
	}
	PanMean = (double)sumvalue[nRastercount]/(width*height);
	for (int k = 0;k < nRastercount;k++)
	{
		MSMean[k]=(double)sumvalue[k]/(width*height);
	}
	delete []sumvalue;
	sumvalue = NULL;
	//2.计算相关系数(全色与多光谱各个波段的相关系数)
	double *relation = new double[nRastercount];
	memset(relation,0,sizeof(double)*nRastercount);
	double temp1=0,temp2=0,temp3=0;

	for (int k = 0;k < nRastercount;k++)
	{
		for (int i = 0;i < height;i++)
		{
			for (int j = 0;j < width;j++)
			{
				temp1 += (pPANImage[i*width+j] - PanMean)*(pMSImage[k][i*width+j] - MSMean[k]);
				temp2 += (pPANImage[i*width+j] - PanMean)*(pPANImage[i*width+j] - PanMean);
				temp3 += (pMSImage[k][i*width+j] - MSMean[k])*(pMSImage[k][i*width+j] - MSMean[k]);
			}
		}
		relation[k] = temp1/sqrt(temp2*temp3);
		temp1 = 0;
		temp2 = 0;
		temp3 = 0;
	}
	delete []MSMean;
	MSMean = NULL;

	unsigned char **FusionImg = new unsigned char *[nRastercount];
	for (int k = 0;k < nRastercount;k++)
	{
		FusionImg[k] = new unsigned char[width*height];
		memset(FusionImg[k],0,sizeof(unsigned char)*width*height);
	}

	for (int i = 0;i < height;i++)
	{
		for (int j = 0;j < width;j++)
		{
			for (int k = 0;k < nRastercount;k++)
			{
				FusionImg[k][i*width+j] = (unsigned char)(((1+abs(relation[k]))*pPANImage[i*width+j] + (1-abs(relation[k]))*pMSImage[k][i*width+j])/2);
			}
		}
	}
	linearStretch(FusionImg,width,height,nRastercount);
	SaveImage(FusionImg,savefusionpath,width,height,nRastercount);

	for (int k = 0; k < nRastercount;k++)
	{
		delete []FusionImg[k];
		FusionImg[k] = NULL;
	}
	delete []FusionImg;
	FusionImg = NULL;

	return ;
}

//////////////////////////////////////////////////////////////////////////
/*                               主成分分析融合方法                     */
//////////////////////////////////////////////////////////////////////////
void principalFusion(const char * filepath1,const char * filepath2,char* savefusionpath)  
{  
	int width =0;
	int height = 0;
	int nRastercount=0;
	//读入多光谱影像
	unsigned char **pMSImage = GDALreadImage(width,height,nRastercount,filepath1);
	/////////////////////////第一步：多光谱影像之间求特征值和特征向量////////////////////////////////////

	//1.1计算波段均值并计算协方差矩阵  
	double* covMatrix = CovMatrix(filepath1);
	//1.2计算协方差所形成的矩阵的特征值与特征向量  
	double eps = 0.000001;   //控制精度要求  
	double *eigenVector = new double[nRastercount * nRastercount];  
	eejcb( covMatrix, nRastercount, eigenVector, eps, 100000);//雅可比法求解  
	//1.3按特征值由大到小的顺序排列特征向量  
	sortEigenVector(eigenVector, covMatrix, nRastercount);
	delete []covMatrix;covMatrix = NULL;
	///////////////////////////第二步：主成分正变换/////////////////////////////////////
	//2.1将多光谱图像进行主分量变换  
	transMatrix(eigenVector,nRastercount,nRastercount);//转置
	double *dPixelX = new double[nRastercount];
	memset(dPixelX,0,sizeof(double)*nRastercount);
	double *dPixelY = new double[nRastercount];
	memset(dPixelY,0,sizeof(double)*nRastercount);
	double ** result = new double*[nRastercount];  //中间计算过程均采用double数据类型，不能用其他类型（其他类型会截断数值，改变结果）
	for ( int i = 0; i < nRastercount; i++ )  
	{  
		result[i] = new double[width*height];  
		for ( int j = 0; j < width*height; j++)  
		{  
			result[i][j] = 0;  //初始化
		}  
	}   
	unsigned char* pHRImage = OpenTif(filepath2);

	//2.2将高光谱分辨率的图像进行主分量变换
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			for(int k=0;k<nRastercount;k++)
			{
				dPixelX[k] = pMSImage[k][i * width + j];
			}
			//主分量变换
			multimatrix(eigenVector,dPixelX,nRastercount,nRastercount,1,dPixelY);//矩阵相乘
			dPixelY[0] = pHRImage[i*width+j];//将高分辨率图像替代第一主分量
			for(int k=0;k<nRastercount;k++)
			{
				result[k][i*width+j] = dPixelY[k];
			}
		}
	}
	int length = ColorLen - 1;
	double *dMax = new double[nRastercount];
	double *dMin = new double[nRastercount];
	for (int i = 0;i < nRastercount;i++)
	{
		dMax[i] = -length;
		dMin[i] = length;
	}
	//统计极值
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			for(int k=0;k<nRastercount;k++)
			{
				double dImageBits = (double)result[k][i*width+j];
				if(dMax[k]<(dImageBits))
					dMax[k] = dImageBits;
				if(dMin[k]>(dImageBits))
					dMin[k] = dImageBits;
			}
		}
	}
	//线性拉伸
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			for(int k=0;k<nRastercount;k++)
			{
				double dImageBits = (double)result[k][i*width+j];
				if((dMax[k]-dMin[k])<=length)
					result[k][i*width+j] = (dImageBits-dMin[k]);
				else
					result[k][i*width+j] = (dImageBits-dMin[k])*length/(dMax[k]-dMin[k]);
			}			
		}
	}

	////////////////////////第三步：进行主分量逆变换//////////////////////////	
	inverseMatrix(eigenVector,nRastercount);//求特征向量矩阵的逆矩阵。	
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			for(int k=0;k<nRastercount;k++)
			{
				dPixelY[k] = result[k][i * width + j];
			}		
			multimatrix(eigenVector,dPixelY,nRastercount,nRastercount,1,dPixelX);
			for(int k=0;k<nRastercount;k++)
			{
				result[k][i*width+j] = dPixelX[k];//数据类型问题将会导致结果影像异常，直接用结果result1来接dPixelX转成unsigned char将有色彩溢出？？
			}//例如PCAerror1
		}
	}
	delete []eigenVector;
	delete []dPixelX;
	delete []dPixelY;
	
	for (int i = 0;i < nRastercount;i++)
	{
		dMax[i] = -length;
		dMin[i] = length;
	}
	//统计极值
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			for(int k=0;k<nRastercount;k++)
			{
				double dImageBits = (double)result[k][i*width+j];
				if(dMax[k]<(dImageBits))
					dMax[k] = dImageBits;
				if(dMin[k]>(dImageBits))
					dMin[k] = dImageBits;
			}
		}
	}
	//线性拉伸
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			for(int k=0;k<nRastercount;k++)
			{
				double dImageBits = (double)result[k][i*width+j];
				if((dMax[k]-dMin[k])<=length)
					result[k][i*width+j] = (dImageBits-dMin[k]);
				else
					result[k][i*width+j] = (dImageBits-dMin[k])*length/(dMax[k]-dMin[k]);
			}			
		}
	}
	delete []dMax;
	delete []dMin;
	unsigned char** result1 = new unsigned char*[nRastercount];//写出影像，采用unsigned char类型  
	for ( int i = 0; i < nRastercount; i++ )  
	{  
		result1[i] = new unsigned char[width*height];  
		for ( int j = 0; j < width*height; j++ )  
		{  
			//result1[i][j] = (unsigned char)result[nRastercount-i-1][j];//???这里为什么需要颠倒波段顺序//顺序颠倒？？？？具体原因为影像数据读入是倒着读的
			result1[i][j] = (unsigned char)result[i][j];
		}  //不颠倒顺序产生的影像色彩不一样//例如PCAerror2
	} 
	//主成分分析存在的问题：
	//倒着读倒着写--PCA1---？
	//顺序读顺序写--PCA2---？
	//顺序读倒着写--PCA3---错
	//倒着读顺着写--PCA4---错
	//为什么会出现这样的情况？
	for (int i = 0;i < nRastercount;i++)
	{
		delete []result[i];
		result[i] = NULL;
	}
	delete []result;
	result = NULL;
	//linearStretch(result1,width,height,nRastercount);
	SaveImage(result1,savefusionpath,width,height,nRastercount);
	for (int i = 0;i < nRastercount;i++)
	{
		delete []result1[i];
		result1[i] = NULL;
	}
	delete []result1;
	result1 = NULL;
    return ;
}  


//////////////////////////////////////////////////////////////////////////
/*                       PrincipalFusion                                */
//////////////////////////////////////////////////////////////////////////
/// <summary>  
/// 矩阵转置.  
/// </summary>  
/// <param name="matrix">矩阵数组指针.</param>  
/// <param name="m">矩阵行数.</param>  
/// <param name="n">矩阵列数.</param>  
void transMatrix( double *matrix, int m, int n )  
{  
	double*p = new double[m * n];  
	if( p == NULL) 
	{
		return;  
	}  
	for( int i = 0; i < m; i++ )  
	{
		for( int j = 0; j < n; j++ )  
		{  
			*( p + i * n + j ) = *( matrix + j * m + i );  //*p[i*n+j]=*matrix[j*m+i];
		} 
	}
	for( int i = 0; i < m; i++ ) 
	{
		for( int j = 0; j < n; j++ )  
		{  
			*( matrix + i * n + j ) = *( p + i * n + j );  
		} 
	}	
	delete []p;  
}  
//矩阵相乘
void multimatrix(double *a, double *b, int m, int n, int k, double *c)
{
	int i,j,l,u;
	for (i=0; i<=m-1; i++)
	{
		for (j=0; j<=k-1; j++)
		{ 
			u=i*k+j; 
			c[u]=0.0;
			for (l=0; l<=n-1; l++)
			{
				c[u]=c[u]+a[i*n+l]*b[l*k+j];
			}
		}
	}
	return;
}
/// <summary>  
/// 矩阵求逆.  
/// </summary>  
/// <param name="matrix">矩阵数组指针.</param>  
/// <param name="n">矩阵阶数.</param>  
void inverseMatrix( double *matrix, int n )  
{  
	int *is, *js, i, j, k, l, u, v;  
	double d, p;  
	is = new int[n * sizeof( int )];  
	js = new int[n * sizeof( int )];  
	for ( k = 0; k <= n - 1; k++ )  
	{  
		d = 0.0;  
		for ( i = k; i <= n - 1; i++ ) 
		{
			for ( j = k; j <= n - 1; j++ )  
			{  
				l = i * n + j;  
				p = fabs( matrix[l] );  
				if ( p > d ) 
				{ 
					d = p; 
					is[k] = i; 
					js[k] = j;
				}
			}
		}		
		if ( d + 1.0 == 1.0 )  //判断该矩阵是否满秩
		{ 
			printf( "err**not inv矩阵求逆失败\n" );
			delete []is;  
			delete []js;  
			return;  
		}  
		if ( is[k] != k )  
		{
			for ( j = 0; j <= n - 1; j++ )  
			{  
				u = k * n + j;  
				v = is[k] * n + j;  
				p = matrix[u];  
				matrix[u] = matrix[v];  
				matrix[v] = p;  
			}
		}
		if ( js[k] != k )  
		{
			for ( i = 0; i <= n - 1; i++ )  
			{  
				u = i * n + k;  
				v = i * n + js[k];  
				p = matrix[u];  
				matrix[u] = matrix[v];  
				matrix[v] = p;  
			} 
		}
		l = k * n + k;  
		matrix[l] = 1.0 / matrix[l];  
		for ( j = 0; j <= n - 1; j++ ) 
		{
			if ( j != k )  
			{ 
				u = k * n + j;
				matrix[u] = matrix[u] * matrix[l];
			}
		}
		for ( i = 0; i <= n - 1; i++ )  
		{
			if ( i != k ) 
			{
				for ( j = 0; j <= n - 1; j++ )
				{
					if ( j != k )  
					{  
						u = i * n + j;  
						matrix[u] = matrix[u] - matrix[i * n + k] * matrix[k * n + j];  
					} 
				}
			}
		}
		for ( i = 0; i <= n - 1; i++ ) 
		{
			if ( i != k )  
			{ 
				u = i * n + k;
				matrix[u] = -matrix[u] * matrix[l];
			}  
		}

	}  
	for ( k = n - 1; k >= 0; k-- )  
	{  
		if ( js[k] != k ) 
		{
			for ( j = 0; j <= n - 1; j++ )  
			{  
				u = k * n + j;  
				v = js[k] * n + j;  
				p = matrix[u];  
				matrix[u] = matrix[v];  
				matrix[v] = p;  
			}  
		}
		if ( is[k] != k )  
		{
			for ( i = 0; i <= n - 1; i++ )  
			{  
				u = i * n + k;  
				v = i * n + is[k];  
				p = matrix[u];  
				matrix[u] = matrix[v];  
				matrix[v] = p;  
			}  
		}
	}  
	delete []is;  
	delete []js;  
}  
/// <summary>  
/// 利用雅格比(Jacobi)方法求实对称矩阵的全部特征值及特征向量.  
/// </summary>  
/// <param name="a">长度为n*n的数组，存放实对称矩阵，返回时对角线存放n个特征值.</param>  
/// <param name="n">矩阵的阶数.</param>  
/// <param name="v">长度为n*n的数组，返回特征向量(按列存储).</param>  
/// <param name="eps">控制精度要求.</param>  
/// <param name="jt">整型变量，控制最大迭代次数.</param>  
/// <returns>返回false表示超过迭代jt次仍未达到精度要求，返回true表示正常返回.</returns>  
bool eejcb( double a[], int n, double v[], double eps, int jt ) 
{  
    int i, j, p, q, u, w, t, s, l;  
    double fm, cn, sn, omega, x, y, d;  
      
    l = 0;  
    //初始化特征向量矩阵使其全为0  
	//相当于对角元素aij（i=j）均为1，而其他元素aij（i≠j）均为0。
    for( i = 0; i <= n - 1; i++ )  
    {  
        v[i * n + i] = 1.0; 
		/************************************************************************/
		/*
		将n*n的矩阵线性化，数组存储均是线性存储，那么矩阵中元素标号为：0,1,2,3,···,n*n-1
		则，对角元素所在位置为;i*n+i,其他的元素位置则表示为：i*n+j(i≠j).
		下面对于多维数组而言。。不能乱用，要按参数设置的来构造！！！
		//*(*(v+i)+i)=1.0; //第i行第i列的元素初始化为1.0
		 //*(*(v+i)+j)=0.0;//第i行第j列的元素初始化为0.0 
		*/                                                                    
		/************************************************************************/
        for( j = 0; j <= n - 1; j++ )  
        {  
            if( i != j ) 
			{
                v[i * n + j] = 0.0; 
			}
			//printf("%f \n",v[i*n+j]);//查看初始化后特征向量矩阵各元素
        }  	
    }  

    while( true ) //循环  
    {  
        fm = 0.0;  
        for( i = 0; i <= n - 1; i++ )   // 求出, 矩阵a( 特征值 ), 中除对角线外其他元素的最大绝对值  
        {  
            //这个最大值是位于a[p][q] ,等于fm  
            for( j = 0; j <= n - 1; j++ )  
            {  
                d = fabs( a[i * n + j] );                  
                if( ( i != j ) && ( d > fm ) )  
                {  
                    fm = d;  
                    p = i;  
                    q = j;  
                }  
            }  
        }  
          
        if( fm < eps )   //精度复合要求  
            return true; //正常返回  
              
        if( l > jt )     //迭代次数太多  
            return false;//失败返回  
              
        l ++;       //   迭代计数器 
		//cout<<"求特征向量矩阵迭代次数："<<l<<endl;
        u = p * n + q;  
        w = p * n + p;  
        t = q * n + p;  
        s = q * n + q;  
        x = -a[u];  
        y = ( a[s] - a[w] ) / 2.0;      //x y的求法不同  
        omega = x / sqrt( x * x + y * y );  //sin2θ  
          
        //tan2θ=x/y = -2.0*a[u]/(a[s]-a[w])  
        if( y < 0.0 )  
            omega = -omega;  
              
        sn = 1.0 + sqrt( 1.0 - omega * omega );  
        sn = omega / sqrt( 2.0 * sn );      //sinθ  
        cn = sqrt( 1.0 - sn * sn );         //cosθ  
          
        fm = a[w];   //   变换前的a[w]   a[p][p]  
        a[w] = fm * cn * cn + a[s] * sn * sn + a[u] * omega;  
        a[s] = fm * sn * sn + a[s] * cn * cn - a[u] * omega;  
        a[u] = 0.0;  
        a[t] = 0.0;  
          
        //   以下是旋转矩阵,旋转了了p行,q行,p列,q列  
        //   但是四个特殊点没有旋转(这四个点在上述语句中发生了变化)  
        //   其他不在这些行和列的点也没变  
        //   旋转矩阵,旋转p行和q行  
        for( j = 0; j <= n - 1; j++ )  
        {  
            if( ( j != p ) && ( j != q ) )  
            {  
                u = p * n + j;  
                w = q * n + j;  
                fm = a[u];  
                a[u] = a[w] * sn + fm * cn;  
                a[w] = a[w] * cn - fm * sn;  
            }  
        }  
          
        //旋转矩阵,旋转p列和q列  
        for( i = 0; i <= n - 1; i++ )  
        {  
            if( ( i != p ) && ( i != q ) )  
            {  
                u = i * n + p;  
                w = i * n + q;  
                fm = a[u];  
                a[u] = a[w] * sn + fm * cn;  
                a[w] = a[w] * cn - fm * sn;  
            }  
        }  
          
        //记录旋转矩阵特征向量  
        for( i = 0; i <= n - 1; i++ )  
        {  
            u = i * n + p;  
            w = i * n + q;  
            fm = v[u];  
            v[u] = v[w] * sn + fm * cn;  
            v[w] = v[w] * cn - fm * sn;  
        } 
    }  
    return true;  
}  

/// <summary>  
/// 线性拉伸.  
/// </summary>  
/// <param name="pResult">图像矩阵.</param>  
/// <param name="width">图像宽.</param>  
/// <param name="height">图像高.</param>  
/// <param name="bandCount">图像波段数.</param>  
void linearStretch( unsigned char** pResult, int width, int height, int bandCount )  
{   
	double thre = ColorLen - 1;
	//对于8bit数据一个像元灰度级最多为256，对于16bit的则是65536
	for ( int i = 0; i < bandCount; i++)  
	{  
		double dMaxValue = -thre, dMinValue = thre;  
		for ( int index = 0; index < width * height; index++ )  
		{  
			if ( dMaxValue < pResult[i][index] )  
			{  
				dMaxValue = pResult[i][index];  
			}  
			if ( dMinValue > pResult[i][index])  
			{  
				dMinValue = pResult[i][index];  
			}  
		}  
		for( int j = 0; j < width * height; j++ )  
		{  
			if ( dMaxValue - dMinValue <= thre)  
			{  
				pResult[i][j] = pResult[i][j] - dMinValue;  
			}  
			else  
			{  
				pResult[i][j] = ( pResult[i][j] - dMinValue ) * thre / (dMaxValue - dMinValue);  
			}  
		}  
	} 
}  
void linearStretch( double** pResult, int width, int height, int bandCount )  
{   
	double thre = ColorLen - 1;
	//对于8bit数据一个像元灰度级最多为256，对于16bit的则是65536
	for ( int i = 0; i < bandCount; i++)  
	{  
		double dMaxValue = -thre, dMinValue = thre;  
		for ( int index = 0; index < width * height; index++ )  
		{  
			if ( dMaxValue < pResult[i][index] )  
			{  
				dMaxValue = pResult[i][index];  
			}  
			if ( dMinValue > pResult[i][index])  
			{  
				dMinValue = pResult[i][index];  
			}  
		}  
		for( int j = 0; j < width * height; j++ )  
		{  
			if ( dMaxValue - dMinValue <= thre)  
			{  
				pResult[i][j] = pResult[i][j] - dMinValue;  
			}  
			else  
			{  
				pResult[i][j] = ( pResult[i][j] - dMinValue ) * thre / (dMaxValue - dMinValue);  
			}  
		}  
	} 
} 

double* CovMatrix(const char* filepath)
{
	//注册、读取图像
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");//使之支持中文路径
	GDALDataset *poDataset = NULL;
	poDataset = (GDALDataset*)GDALOpen(filepath,GA_ReadOnly);
	if(poDataset == NULL)
	{
		cout<<"无法打开图像！"<<endl;
		GDALDestroyDriverManager();
	}
	//获取图像数据的参数
	int width = poDataset->GetRasterXSize();
	int height = poDataset->GetRasterYSize();
	int nRastercount = poDataset->GetRasterCount();
	//计算多光谱图像波段均值
	double* bandMean = new double [nRastercount];  
	for ( int i = 0; i < nRastercount; i++ )  
	{  
		double dMaxValue, dMinValue;  
		poDataset->GetRasterBand( nRastercount - i)->ComputeStatistics( FALSE, &dMinValue, &dMaxValue, bandMean + i, 0, NULL, NULL );  
	}  
	if ( bandMean == NULL )  
	{  
		cout << "统计波段均值失败！" << endl;  
		return NULL;  
	} 

	//求图像协方差矩阵
	double *Covariance = new double[nRastercount * nRastercount];  
	memset(Covariance,0,sizeof(double)*nRastercount*nRastercount);
	int index = 0;  
	for ( int i = 0; i < nRastercount; i++ )  
	{  
		unsigned char* poData1 = new unsigned char[ height * width];  
		int bandList = {nRastercount - i};  
		poDataset->RasterIO( GF_Read, 0, 0, width, height, poData1, width, height, GDT_Byte, 1, &bandList, 0, 0, 0 );  

		for ( int j = 0; j < nRastercount; j++ )  
		{  
			unsigned char* poData2 = new unsigned char[ height * width];  
			int bandList = {nRastercount - j};  
			poDataset->RasterIO( GF_Read, 0, 0, width, height, poData2, width, height, GDT_Byte, 1, &bandList, 0, 0, 0  );  
			double sum = 0;  

			for ( int pix = 0; pix < height * width; pix++ )  
			{  
				sum += ( poData1[pix] - bandMean[i] ) * ( poData2[pix] - bandMean[j] );  
			}  
			Covariance[index++] = sum * 1.0 / ( height * width - 1 );  
			delete []poData2;poData2 = NULL;
		}  
		delete []poData1;poData1 = NULL;
	}
	GDALClose(poDataset);
	poDataset = NULL;
	delete []bandMean;bandMean = NULL;
	return Covariance;
}

/// <summary>  
/// 按特征值大小排列特征向量.  
/// </summary>  
/// <param name="eigenVector">特征向量.</param>  
/// <param name="covAfterEejcb">经过Jacobi方法返回的原协方差矩阵，对角线为特征值.</param>  
void sortEigenVector( double * eigenVector, double * covAfterEejcb,int nRastercount)  
{  
	for( int i = 0; i < nRastercount - 1; i++ )  
	{  
		for( int j = i + 1; j < nRastercount; j++ )  
		{  
			if( covAfterEejcb[j * nRastercount + j] > covAfterEejcb[i * nRastercount + i] )  
			{  
				double temp = 0;  
				temp = covAfterEejcb[j * nRastercount + j];  
				covAfterEejcb[j * nRastercount + j] = covAfterEejcb[i * nRastercount + i];  
				covAfterEejcb[i * nRastercount + i] = temp;  
				for( int k = 0; k < nRastercount; k++ )  
				{  
					temp = 0;  
					temp = eigenVector[k * nRastercount + j];  
					eigenVector[k * nRastercount + j] = eigenVector[k * nRastercount + i];  
					eigenVector[k * nRastercount + i] = temp;  
				}  
			}  
		}  
	}  
} 
/*
/// <summary>  
/// PCA变换.  
/// </summary>  
/// <param name="imgMatrix">图像矩阵.</param>  
unsigned char** PCATransform( unsigned char **imgMatrix, double  *eigenVector,int width,int height,int nRastercount)  
{  
	double **vec = new double*[nRastercount];  
	for ( int i = 0; i < nRastercount; i++ )  
	{  
		vec[i] = new double[nRastercount];  
		for ( int j = 0; j < nRastercount; j++)  
		{  
			vec[i][j] = eigenVector[i + j * nRastercount];  
		}  
	}  //写成列向量

	// 构造结果矩阵  
	unsigned char** result = new unsigned char*[nRastercount];  
	for ( int i = 0; i < nRastercount; i++ )  
	{  
		result[i] = new unsigned char[width*height];  
		for ( int j = 0; j < width*height ; j++ )  
		{  
			result[i][j] = 0;  //初始化
		}  
	}  

	// 矩阵相乘  
	for ( int i = 0; i < nRastercount; i++ )  
	{  
		for ( int j = 0; j < width*height; j++ )  
		{  
			for ( int k = 0; k < nRastercount; k++ )  
			{  
				result[i][j] +=  imgMatrix[k][j] * vec[i][k];  
			}  
		}  
	} 
	for (int i = 0; i < nRastercount;i++)
	{
		delete []vec[i];
		vec[i] = NULL;
	}
	delete []vec;
	return result;  
}  
/// <summary>  
/// PCA逆变换.  
/// </summary>  
/// <param name="imgMatrix">图像矩阵.</param>  
unsigned char** inversePCA( unsigned char **imgMatrix, double *eigenVector,int width,int height,int n)  
{  
	transMatrix(eigenVector,n,n);
	// 求特征向量矩阵的逆矩阵  
	inverseMatrix(eigenVector, n);  
	for (int i=0;i<9;i++)
	{
		printf("特征向量求逆：%.4f\n",eigenVector[i]);
	}
	unsigned char** resAfterInversPCA = PCATransform( imgMatrix, eigenVector,width,height,n);  
	for (int i = 0; i < sqrt((double)n);i++)
	{
		delete []imgMatrix[i];
		imgMatrix[i] = NULL;
	}
	imgMatrix = NULL;
	return resAfterInversPCA;  
} 

 //<summary>  
 //PCA融合.  
 //</summary>  
unsigned char**principalFusion(const char * filepath1,const char * filepath2)  
{  
	int width =0;
	int height = 0;
	int nRastercount=0;
	//读入多光谱影像
	unsigned char **pMSImage = GDALreadImage(width,height,nRastercount,filepath1);
	/////////////////////////第一步：多光谱影像之间求特征值和特征向量////////////////////////////////////
	//int width*height = width*height;
    //1.1计算波段均值并计算协方差矩阵  
    double* covMatrix = CovMatrix(filepath1);
    //1.2计算协方差所形成的矩阵的特征值与特征向量  
    double eps = 0.000001;   //控制精度要求  
    double *eigenVector = new double[nRastercount * nRastercount];  
    eejcb( covMatrix, nRastercount, eigenVector, eps, 100000);//雅可比法求解  
    //1.3按特征值由大到小的顺序排列特征向量  
    sortEigenVector(eigenVector, covMatrix, nRastercount);
	delete []covMatrix;covMatrix = NULL;
	///////////////////////////第二步：主成分正变换/////////////////////////////////////
    //2.1将多光谱图像进行主分量变换  
	transMatrix(eigenVector,nRastercount,nRastercount);//转置
	double *dPixelX = new double[nRastercount];
	memset(dPixelX,0,sizeof(double)*nRastercount);
	double *dPixelY = new double[nRastercount];
	memset(dPixelY,0,sizeof(double)*nRastercount);
	double ** result = new double*[nRastercount];  //中间计算过程均采用double数据类型，不能用其他类型（其他类型会截断数值，改变结果）
	for ( int i = 0; i < nRastercount; i++ )  
	{  
		result[i] = new double[width*height];  
		for ( int j = 0; j < width*height ; j++)  
		{  
			result[i][j] = 0;  //初始化
		}  
	}   
	unsigned char* pHRImage = OpenTif(filepath2);

	//2.2将高光谱分辨率的图像进行主分量变换
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			for(int k=0;k<nRastercount;k++)
			{
				dPixelX[k] = pMSImage[k][i * width + j];
			}
			//主分量变换
			multimatrix(eigenVector,dPixelX,nRastercount,nRastercount,1,dPixelY);//矩阵相乘
			dPixelY[0] = pHRImage[i*width+j];//将高分辨率图像替代第一主分量
			for(int k=0;k<nRastercount;k++)
			{
				result[k][i*width+j] = dPixelY[k];
			}
		}
	}
	int length = ColorLen - 1;
	double *dMax = new double[nRastercount];
	double *dMin = new double[nRastercount];
	for (int i = 0;i < nRastercount;i++)
	{
		dMax[i] = -length;
		dMin[i] = length;
	}
	//统计极值
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			for(int k=0;k<nRastercount;k++)
			{
				double dImageBits = (double)result[k][i*width+j];
				if(dMax[k]<(dImageBits))
					dMax[k] = dImageBits;
				if(dMin[k]>(dImageBits))
					dMin[k] = dImageBits;
			}
		}
	}
	//线性拉伸
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			for(int k=0;k<nRastercount;k++)
			{
				double dImageBits = (double)result[k][i*width+j];
				if((dMax[k]-dMin[k])<length)
					result[k][i*width+j] = (dImageBits-dMin[k]);
				else
					result[k][i*width+j] = (dImageBits-dMin[k])*length/(dMax[k]-dMin[k]);
			}			
		}
	}
	////////////////////////第三步：进行主分量逆变换//////////////////////////
	inverseMatrix(eigenVector,nRastercount);//求特征向量矩阵的逆矩阵。	
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			for(int k=0;k<nRastercount;k++)
			{
				dPixelX[k] = result[k][i * width + j];
			}		
			multimatrix(eigenVector,dPixelX,nRastercount,nRastercount,1,dPixelY);
			for(int k=0;k<nRastercount;k++)
			{
				result[k][i*width+j] = dPixelY[k];
			}
		}
	}
	delete []eigenVector;
	delete []dPixelX;
	delete []dPixelY;
	//重新赋值
	for (int i = 0;i < nRastercount;i++)
	{
		dMax[i] = -length;
		dMin[i] = length;
	}
	//统计极值
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			for(int k=0;k<nRastercount;k++)
			{
				double dImageBits = (double)result[k][i*width+j];
				if(dMax[k]<(dImageBits))
					dMax[k] = dImageBits;
				if(dMin[k]>(dImageBits))
					dMin[k] = dImageBits;
			}
		}
	}
	//线性拉伸
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			for(int k=0;k<nRastercount;k++)
			{
				double dImageBits = (double)result[k][i*width+j];
				if((dMax[k]-dMin[k])<length)
					result[k][i*width+j] = (dImageBits-dMin[k]);
				else
					result[k][i*width+j] = (dImageBits-dMin[k])*length/(dMax[k]-dMin[k]);
			}			
		}
	}
	delete []dMax;
	delete []dMin;
	unsigned char** result1 = new unsigned char*[nRastercount];//写出影像，采用unsigned char类型  
	for ( int i = 0; i < nRastercount; i++ )  
	{  
		result1[i] = new unsigned char[width*height];  
		for ( int j = 0; j < width*height ; j++ )  
		{  
			//result1[i][j] = (unsigned char)result[nRastercount-i-1][j];//???这里为什么需要颠倒波段顺序
			result1[i][j] = (unsigned char)result[i][j];
		}  
	} 

	for (int i = 0;i < nRastercount;i++)
	{
		delete []result[i];
		result[i] = NULL;
	}
	delete []result;
	result = NULL;
	return result1;	
}  


*/