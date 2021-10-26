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
/*                              HPF����ͨ�˲���                         */
//////////////////////////////////////////////////////////////////////////
void HPFFusion(const char* msfilepath,const char* panfilepath,char* savefusionpath)  
{  
	int width =0;
	int height = 0;
	int nRastercount=0;
	//��������Ӱ��
	unsigned char **pMSImage = GDALreadImage(width,height,nRastercount,msfilepath);
	//�˲�����ѡ��3��3
	int FilterWidth = 3;
	//����ȫɫӰ��
	unsigned char *pPANImage = OpenTif(panfilepath);
	//�м����
	double *temppan = new double[width*height];
	memset(temppan,0,sizeof(double)*width*height);
	//�˲���
	unsigned char *temp = new unsigned char[FilterWidth*FilterWidth];
	memset(temp,0,sizeof(unsigned char)*FilterWidth*FilterWidth);

	//i,j��1��ʼ��width-1����������i*width+j���ܱ�ʾ�˲����м�ֵ���൱����ΧһȦ������ԭֵ���ڲ������˲����¸�ֵ
	for (int i = FilterWidth/2;i < height - FilterWidth/2;++i)       
	{
		for (int j = FilterWidth/2;j < width - FilterWidth/2;++j)
		{
			int index=0;double sum = 0;
			for (int k = -FilterWidth/2;k <= FilterWidth/2;++k)
			{
				for (int h = -FilterWidth/2;h <= FilterWidth/2;++h)
				{
					temp[index] = pPANImage[(i + k) * width + (j + h)];//��������˲����и����Ҷ�ֵ
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
	//�Ѹ�Ƶ�ɷּӵ������ͼ������
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
			//ֻ��ĳһ�������ϼӸ�Ƶ��Ϣ�����Ӱ����ģ���ģ����μ�ֱ��ʲ�һ�£������ܴﵽ��ǿӰ��ֱ��ʵ�Ҫ��
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
	//ͬʱ�Ӹ�Ƶ��Ϣ���������죬�����ε�����ֵ������ƽ���Ӷ�ʹ���ں�Ӱ��׻�����ʹ������Ϣ��ʧ����
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
/*                    �˻��Ա任�ں��㷨��Multiply�Լ��Ľ��㷨��        */
//////////////////////////////////////////////////////////////////////////
void MultiplyFusion(const char* msfilepath,const char* panfilepath,char* savefusionpath,int flag)
{
	int width =0;
	int height = 0;
	int nRastercount=0;
	double RGBsum = 0.0;
	double panvalue = 0.0;
	//��������Ӱ��
	unsigned char **pMSImage = GDALreadImage(width,height,nRastercount,msfilepath);
	//����ȫɫӰ��
	unsigned char* pPANImage = OpenTif(panfilepath);
	//�м����
	double**tempimg = new double*[nRastercount];
	for (int k = 0;k < nRastercount;k++)
	{
		tempimg[k] = new double[width*height];
		memset(tempimg[k],0,sizeof(double)*width*height);
	}
	switch(flag)
	{
	case 0:
		cout<<"�˻��Ա任"<<endl;	
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
		cout<<"���Ը�����˻��任"<<endl;
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
/*                   ��ɫ�任�ں��㷨��IHS�ںϷ����Լ�����㷨��        */
//////////////////////////////////////////////////////////////////////////
void IHSFusion(const char* msfilepath,const char* panfilepath,char* savefusionpath,int flag)
{
	int width =0;
	int height = 0;
	int nRastercount=0;
	//��������Ӱ��
	unsigned char **pMSImage = GDALreadImage(width,height,nRastercount,msfilepath);
	if (nRastercount != 3)
	{
		cout<<"����ĵͷֱ���Ӱ�񲨶������󣬲��ܽ��и��ںϲ�����"<<endl;
		for (int k = 0; k < nRastercount;k++)
		{
			delete []pMSImage[k];
			pMSImage[k] = NULL;
		}
		delete []pMSImage;
		pMSImage = NULL;
		return;
	}
	//����ȫɫӰ��
	unsigned char* pPANImage = OpenTif(panfilepath);
	//�м����
	double**tempimg = new double*[nRastercount];//�����м���̣����Ӱ���Ч���仯���󣿣���
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
		cout<<"IHS�任"<<endl;	
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
				//�������任
				fI = RGBsum/3;
				max_min_value(dRGB,nRastercount,dmax,dmin);
				fS = (RGBsum - 3*dmin)/RGBsum;
				double angle = acos(0.5*(2*dRGB[0]-dRGB[1]-dRGB[2])/(sqrt((dRGB[0]-dRGB[1])*(dRGB[0]-dRGB[1])+(dRGB[0]-dRGB[2])*(dRGB[1]-dRGB[2]))));
				angle = angle*180.0/PI;
				if (dRGB[2]>dRGB[1])//BG�Ƚ�
				{
					fH = (360-angle)/360;
				}
				else
				{
					fH = angle/360;
				}
				
				dRGBt[0] = fI;dRGBt[1] = fS;dRGBt[2] = fH;
				//�ø߿ռ�ֱ��ʵ�ͼ���滻Intensity���ȷ���
				dRGBt[0] = (double)pPANImage[i*width+j];
				//������任
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
				//�����ڴ���
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
		cout<<"��������IHS�任"<<endl;	
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
				//�������任
				fI = RGBsum/3;
				max_min_value(dRGB,nRastercount,dmax,dmin);
				fS = (RGBsum - 3*dmin)/RGBsum;
				if (dmin == dRGB[0])//R��С
				{
					//fS = (RGBsum - 3*dRGB[0])/RGBsum;
					if((RGBsum-3*dRGB[0]) != 0) 
					{
						fH = (dRGB[2]-dRGB[0])/(RGBsum-3*dRGB[0])+1;
					}
					else //���ж�Ԥ��������������һ�㶼��true
					{
						fH = (dRGB[2]-dRGB[0])/(RGBsum-3*dRGB[0]+1)+1;	
					}
				}
				if (dmin == dRGB[2])//B��С
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
				if (dmin == dRGB[1])//G��С
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
				//�ø߿ռ�ֱ��ʵ�ͼ���滻Intensity���ȷ���
				dRGBt[0] = (double)pPANImage[i*width+j];
				//������任
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
				//�����ڴ���
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
		cout<<"������ͨԲ����IHS�任(���IHS)"<<endl;	
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
				//�������任//G��B�Ƚ�
				if(dRGB[1]>=dRGB[2]) C = 0;
				if(dRGB[1]<dRGB[2]) C = PI;
				fI = (1/sqrt(3.0))*(dRGB[2]+dRGB[1]+dRGB[0]);
				v1 = (dRGB[2]+dRGB[1]-2*dRGB[0])/(sqrt(6.0));
				v2 = (dRGB[2]-dRGB[1])/(sqrt(2.0));
				fH = atan(v1/v2) + C;
				fS = sqrt(v1*v1+v2*v2);			   
				dRGBt[0] = fI;dRGBt[1] = fS;dRGBt[2] = fH;
				//�ø߿ռ�ֱ��ʵ�ͼ���滻Intensity���ȷ���
				dRGBt[0] = pPANImage[i*width+j];
				//������任
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
				//�����ڴ���
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
		cout<<"���������Բ����IHS�任"<<endl;
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
				//�������任//G��B�Ƚ�
				if(dRGB[1]>=dRGB[2]) C = 0;
				if(dRGB[1]<dRGB[2]) C = PI;//��
				fI = (1/sqrt(3.0))*(dRGB[2]+dRGB[1]+dRGB[0]);
				v1 = (dRGB[2]+dRGB[1]-2*dRGB[0])/(sqrt(6.0));
				v2 = (dRGB[2]-dRGB[1])/(sqrt(2.0));
				fH = atan(v1/v2) + C;
				fS = sqrt(v1*v1+v2*v2);			   
				dRGBt[0] = fI;dRGBt[1] = fS;dRGBt[2] = fH;
				//�ø߿ռ�ֱ��ʵ�ͼ���滻Intensity���ȷ���
				dRGBt[0] = pPANImage[i*width+j];
				//������任
				br = (1/sqrt(3.0))*dRGBt[0]-(2/sqrt(6.0))*v1;
				bg = (1/sqrt(3.0))*dRGBt[0]+(1/sqrt(6.0))*v1-(1/sqrt(2.0))*v2;
				bb = (1/sqrt(3.0))*dRGBt[0]+(1/sqrt(6.0))*v1+(1/sqrt(2.0))*v2;	
				/*br *= 1.5;
				bg *= 1.5;
				bb *= 1.5;*///????������������
				if(br>255) br = 255;
				if(br<0)   br = 0;
				if(bg>255) bg = 255;
				if(bg<0)  bg = 0;
				if(bb>255) bb = 255;
				if(bb<0)  bb = 0;
				dRGB[0] = br;dRGB[1] = bg;dRGB[2] = bb;
				//�����ڴ���
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

//�����������ֵ����Сֵ
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
/*                ��ֵ�ں��㷨��Brovey�ںϷ����Լ���ظĽ��㷨��        */
//////////////////////////////////////////////////////////////////////////
void BroveyFusion(const char* msfilepath,const char* panfilepath,char* savefusionpath)
{
	int width =0;
	int height = 0;
	int nRastercount=0;
	//��������Ӱ��
	unsigned char **pMSImage = GDALreadImage(width,height,nRastercount,msfilepath);
	//����ȫɫӰ��
	unsigned char* pPANImage = OpenTif(panfilepath);
	//�м����
	double**tempimg = new double*[nRastercount];//�����м���̣����Ӱ���Ч���仯���󣿣���
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
/*                       �������ϵ���ļ�Ȩ�ںϷ���                     */
//////////////////////////////////////////////////////////////////////////
void CorrelationFusion(const char* msfilepath,const char* panfilepath,char* savefusionpath)
{
	int width =0;
	int height = 0;
	int nRastercount=0;
	//��������Ӱ��
	unsigned char **pMSImage = GDALreadImage(width,height,nRastercount,msfilepath);
	//����ȫɫӰ��
	unsigned char* pPANImage = OpenTif(panfilepath);
	//1.��������ξ�ֵ
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
	//2.�������ϵ��(ȫɫ�����׸������ε����ϵ��)
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
/*                               ���ɷַ����ںϷ���                     */
//////////////////////////////////////////////////////////////////////////
void principalFusion(const char * filepath1,const char * filepath2,char* savefusionpath)  
{  
	int width =0;
	int height = 0;
	int nRastercount=0;
	//��������Ӱ��
	unsigned char **pMSImage = GDALreadImage(width,height,nRastercount,filepath1);
	/////////////////////////��һ���������Ӱ��֮��������ֵ����������////////////////////////////////////

	//1.1���㲨�ξ�ֵ������Э�������  
	double* covMatrix = CovMatrix(filepath1);
	//1.2����Э�������γɵľ��������ֵ����������  
	double eps = 0.000001;   //���ƾ���Ҫ��  
	double *eigenVector = new double[nRastercount * nRastercount];  
	eejcb( covMatrix, nRastercount, eigenVector, eps, 100000);//�ſɱȷ����  
	//1.3������ֵ�ɴ�С��˳��������������  
	sortEigenVector(eigenVector, covMatrix, nRastercount);
	delete []covMatrix;covMatrix = NULL;
	///////////////////////////�ڶ��������ɷ����任/////////////////////////////////////
	//2.1�������ͼ������������任  
	transMatrix(eigenVector,nRastercount,nRastercount);//ת��
	double *dPixelX = new double[nRastercount];
	memset(dPixelX,0,sizeof(double)*nRastercount);
	double *dPixelY = new double[nRastercount];
	memset(dPixelY,0,sizeof(double)*nRastercount);
	double ** result = new double*[nRastercount];  //�м������̾�����double�������ͣ��������������ͣ��������ͻ�ض���ֵ���ı�����
	for ( int i = 0; i < nRastercount; i++ )  
	{  
		result[i] = new double[width*height];  
		for ( int j = 0; j < width*height; j++)  
		{  
			result[i][j] = 0;  //��ʼ��
		}  
	}   
	unsigned char* pHRImage = OpenTif(filepath2);

	//2.2���߹��׷ֱ��ʵ�ͼ������������任
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			for(int k=0;k<nRastercount;k++)
			{
				dPixelX[k] = pMSImage[k][i * width + j];
			}
			//�������任
			multimatrix(eigenVector,dPixelX,nRastercount,nRastercount,1,dPixelY);//�������
			dPixelY[0] = pHRImage[i*width+j];//���߷ֱ���ͼ�������һ������
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
	//ͳ�Ƽ�ֵ
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
	//��������
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

	////////////////////////��������������������任//////////////////////////	
	inverseMatrix(eigenVector,nRastercount);//��������������������	
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
				result[k][i*width+j] = dPixelX[k];//�����������⽫�ᵼ�½��Ӱ���쳣��ֱ���ý��result1����dPixelXת��unsigned char����ɫ���������
			}//����PCAerror1
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
	//ͳ�Ƽ�ֵ
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
	//��������
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
	unsigned char** result1 = new unsigned char*[nRastercount];//д��Ӱ�񣬲���unsigned char����  
	for ( int i = 0; i < nRastercount; i++ )  
	{  
		result1[i] = new unsigned char[width*height];  
		for ( int j = 0; j < width*height; j++ )  
		{  
			//result1[i][j] = (unsigned char)result[nRastercount-i-1][j];//???����Ϊʲô��Ҫ�ߵ�����˳��//˳��ߵ�������������ԭ��ΪӰ�����ݶ����ǵ��Ŷ���
			result1[i][j] = (unsigned char)result[i][j];
		}  //���ߵ�˳�������Ӱ��ɫ�ʲ�һ��//����PCAerror2
	} 
	//���ɷַ������ڵ����⣺
	//���Ŷ�����д--PCA1---��
	//˳���˳��д--PCA2---��
	//˳�������д--PCA3---��
	//���Ŷ�˳��д--PCA4---��
	//Ϊʲô����������������
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
/// ����ת��.  
/// </summary>  
/// <param name="matrix">��������ָ��.</param>  
/// <param name="m">��������.</param>  
/// <param name="n">��������.</param>  
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
//�������
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
/// ��������.  
/// </summary>  
/// <param name="matrix">��������ָ��.</param>  
/// <param name="n">�������.</param>  
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
		if ( d + 1.0 == 1.0 )  //�жϸþ����Ƿ�����
		{ 
			printf( "err**not inv��������ʧ��\n" );
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
/// �����Ÿ��(Jacobi)������ʵ�Գƾ����ȫ������ֵ����������.  
/// </summary>  
/// <param name="a">����Ϊn*n�����飬���ʵ�Գƾ��󣬷���ʱ�Խ��ߴ��n������ֵ.</param>  
/// <param name="n">����Ľ���.</param>  
/// <param name="v">����Ϊn*n�����飬������������(���д洢).</param>  
/// <param name="eps">���ƾ���Ҫ��.</param>  
/// <param name="jt">���ͱ�������������������.</param>  
/// <returns>����false��ʾ��������jt����δ�ﵽ����Ҫ�󣬷���true��ʾ��������.</returns>  
bool eejcb( double a[], int n, double v[], double eps, int jt ) 
{  
    int i, j, p, q, u, w, t, s, l;  
    double fm, cn, sn, omega, x, y, d;  
      
    l = 0;  
    //��ʼ��������������ʹ��ȫΪ0  
	//�൱�ڶԽ�Ԫ��aij��i=j����Ϊ1��������Ԫ��aij��i��j����Ϊ0��
    for( i = 0; i <= n - 1; i++ )  
    {  
        v[i * n + i] = 1.0; 
		/************************************************************************/
		/*
		��n*n�ľ������Ի�������洢�������Դ洢����ô������Ԫ�ر��Ϊ��0,1,2,3,������,n*n-1
		�򣬶Խ�Ԫ������λ��Ϊ;i*n+i,������Ԫ��λ�����ʾΪ��i*n+j(i��j).
		������ڶ�ά������ԡ����������ã�Ҫ���������õ������죡����
		//*(*(v+i)+i)=1.0; //��i�е�i�е�Ԫ�س�ʼ��Ϊ1.0
		 //*(*(v+i)+j)=0.0;//��i�е�j�е�Ԫ�س�ʼ��Ϊ0.0 
		*/                                                                    
		/************************************************************************/
        for( j = 0; j <= n - 1; j++ )  
        {  
            if( i != j ) 
			{
                v[i * n + j] = 0.0; 
			}
			//printf("%f \n",v[i*n+j]);//�鿴��ʼ�����������������Ԫ��
        }  	
    }  

    while( true ) //ѭ��  
    {  
        fm = 0.0;  
        for( i = 0; i <= n - 1; i++ )   // ���, ����a( ����ֵ ), �г��Խ���������Ԫ�ص�������ֵ  
        {  
            //������ֵ��λ��a[p][q] ,����fm  
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
          
        if( fm < eps )   //���ȸ���Ҫ��  
            return true; //��������  
              
        if( l > jt )     //��������̫��  
            return false;//ʧ�ܷ���  
              
        l ++;       //   ���������� 
		//cout<<"�����������������������"<<l<<endl;
        u = p * n + q;  
        w = p * n + p;  
        t = q * n + p;  
        s = q * n + q;  
        x = -a[u];  
        y = ( a[s] - a[w] ) / 2.0;      //x y���󷨲�ͬ  
        omega = x / sqrt( x * x + y * y );  //sin2��  
          
        //tan2��=x/y = -2.0*a[u]/(a[s]-a[w])  
        if( y < 0.0 )  
            omega = -omega;  
              
        sn = 1.0 + sqrt( 1.0 - omega * omega );  
        sn = omega / sqrt( 2.0 * sn );      //sin��  
        cn = sqrt( 1.0 - sn * sn );         //cos��  
          
        fm = a[w];   //   �任ǰ��a[w]   a[p][p]  
        a[w] = fm * cn * cn + a[s] * sn * sn + a[u] * omega;  
        a[s] = fm * sn * sn + a[s] * cn * cn - a[u] * omega;  
        a[u] = 0.0;  
        a[t] = 0.0;  
          
        //   ��������ת����,��ת����p��,q��,p��,q��  
        //   �����ĸ������û����ת(���ĸ�������������з����˱仯)  
        //   ����������Щ�к��еĵ�Ҳû��  
        //   ��ת����,��תp�к�q��  
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
          
        //��ת����,��תp�к�q��  
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
          
        //��¼��ת������������  
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
/// ��������.  
/// </summary>  
/// <param name="pResult">ͼ�����.</param>  
/// <param name="width">ͼ���.</param>  
/// <param name="height">ͼ���.</param>  
/// <param name="bandCount">ͼ�񲨶���.</param>  
void linearStretch( unsigned char** pResult, int width, int height, int bandCount )  
{   
	double thre = ColorLen - 1;
	//����8bit����һ����Ԫ�Ҷȼ����Ϊ256������16bit������65536
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
	//����8bit����һ����Ԫ�Ҷȼ����Ϊ256������16bit������65536
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
	//ע�ᡢ��ȡͼ��
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");//ʹ֧֮������·��
	GDALDataset *poDataset = NULL;
	poDataset = (GDALDataset*)GDALOpen(filepath,GA_ReadOnly);
	if(poDataset == NULL)
	{
		cout<<"�޷���ͼ��"<<endl;
		GDALDestroyDriverManager();
	}
	//��ȡͼ�����ݵĲ���
	int width = poDataset->GetRasterXSize();
	int height = poDataset->GetRasterYSize();
	int nRastercount = poDataset->GetRasterCount();
	//��������ͼ�񲨶ξ�ֵ
	double* bandMean = new double [nRastercount];  
	for ( int i = 0; i < nRastercount; i++ )  
	{  
		double dMaxValue, dMinValue;  
		poDataset->GetRasterBand( nRastercount - i)->ComputeStatistics( FALSE, &dMinValue, &dMaxValue, bandMean + i, 0, NULL, NULL );  
	}  
	if ( bandMean == NULL )  
	{  
		cout << "ͳ�Ʋ��ξ�ֵʧ�ܣ�" << endl;  
		return NULL;  
	} 

	//��ͼ��Э�������
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
/// ������ֵ��С������������.  
/// </summary>  
/// <param name="eigenVector">��������.</param>  
/// <param name="covAfterEejcb">����Jacobi�������ص�ԭЭ������󣬶Խ���Ϊ����ֵ.</param>  
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
/// PCA�任.  
/// </summary>  
/// <param name="imgMatrix">ͼ�����.</param>  
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
	}  //д��������

	// ����������  
	unsigned char** result = new unsigned char*[nRastercount];  
	for ( int i = 0; i < nRastercount; i++ )  
	{  
		result[i] = new unsigned char[width*height];  
		for ( int j = 0; j < width*height ; j++ )  
		{  
			result[i][j] = 0;  //��ʼ��
		}  
	}  

	// �������  
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
/// PCA��任.  
/// </summary>  
/// <param name="imgMatrix">ͼ�����.</param>  
unsigned char** inversePCA( unsigned char **imgMatrix, double *eigenVector,int width,int height,int n)  
{  
	transMatrix(eigenVector,n,n);
	// ��������������������  
	inverseMatrix(eigenVector, n);  
	for (int i=0;i<9;i++)
	{
		printf("�����������棺%.4f\n",eigenVector[i]);
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
 //PCA�ں�.  
 //</summary>  
unsigned char**principalFusion(const char * filepath1,const char * filepath2)  
{  
	int width =0;
	int height = 0;
	int nRastercount=0;
	//��������Ӱ��
	unsigned char **pMSImage = GDALreadImage(width,height,nRastercount,filepath1);
	/////////////////////////��һ���������Ӱ��֮��������ֵ����������////////////////////////////////////
	//int width*height = width*height;
    //1.1���㲨�ξ�ֵ������Э�������  
    double* covMatrix = CovMatrix(filepath1);
    //1.2����Э�������γɵľ��������ֵ����������  
    double eps = 0.000001;   //���ƾ���Ҫ��  
    double *eigenVector = new double[nRastercount * nRastercount];  
    eejcb( covMatrix, nRastercount, eigenVector, eps, 100000);//�ſɱȷ����  
    //1.3������ֵ�ɴ�С��˳��������������  
    sortEigenVector(eigenVector, covMatrix, nRastercount);
	delete []covMatrix;covMatrix = NULL;
	///////////////////////////�ڶ��������ɷ����任/////////////////////////////////////
    //2.1�������ͼ������������任  
	transMatrix(eigenVector,nRastercount,nRastercount);//ת��
	double *dPixelX = new double[nRastercount];
	memset(dPixelX,0,sizeof(double)*nRastercount);
	double *dPixelY = new double[nRastercount];
	memset(dPixelY,0,sizeof(double)*nRastercount);
	double ** result = new double*[nRastercount];  //�м������̾�����double�������ͣ��������������ͣ��������ͻ�ض���ֵ���ı�����
	for ( int i = 0; i < nRastercount; i++ )  
	{  
		result[i] = new double[width*height];  
		for ( int j = 0; j < width*height ; j++)  
		{  
			result[i][j] = 0;  //��ʼ��
		}  
	}   
	unsigned char* pHRImage = OpenTif(filepath2);

	//2.2���߹��׷ֱ��ʵ�ͼ������������任
	for(int i=0;i<height;i++)
	{
		for(int j=0;j<width;j++)
		{
			for(int k=0;k<nRastercount;k++)
			{
				dPixelX[k] = pMSImage[k][i * width + j];
			}
			//�������任
			multimatrix(eigenVector,dPixelX,nRastercount,nRastercount,1,dPixelY);//�������
			dPixelY[0] = pHRImage[i*width+j];//���߷ֱ���ͼ�������һ������
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
	//ͳ�Ƽ�ֵ
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
	//��������
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
	////////////////////////��������������������任//////////////////////////
	inverseMatrix(eigenVector,nRastercount);//��������������������	
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
	//���¸�ֵ
	for (int i = 0;i < nRastercount;i++)
	{
		dMax[i] = -length;
		dMin[i] = length;
	}
	//ͳ�Ƽ�ֵ
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
	//��������
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
	unsigned char** result1 = new unsigned char*[nRastercount];//д��Ӱ�񣬲���unsigned char����  
	for ( int i = 0; i < nRastercount; i++ )  
	{  
		result1[i] = new unsigned char[width*height];  
		for ( int j = 0; j < width*height ; j++ )  
		{  
			//result1[i][j] = (unsigned char)result[nRastercount-i-1][j];//???����Ϊʲô��Ҫ�ߵ�����˳��
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