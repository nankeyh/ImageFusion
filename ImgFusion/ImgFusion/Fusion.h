#pragma once

#ifndef   __IMAGEFUSION_H__
#define   __IMAGEFUSION_H__

#include <iostream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include "Function.h"
//#include "Macro.h"


using namespace std;
void HPFFusion(const char* msfilepath,const char* panfilepath,char* savefusionpath) ;
void IHSFusion(const char* msfilepath,const char* panfilepath,char* savefusionpath,int flag);
void MultiplyFusion(const char* msfilepath,const char* panfilepath,char* savefusionpath,int flag);
void BroveyFusion(const char* msfilepath,const char* panfilepath,char* savefusionpath);
void CorrelationFusion(const char* msfilepath,const char* panfilepath,char* savefusionpath);
void principalFusion(const char * filepath1,const char * filepath2,char* savefusionpath);

//////////////////////////////////////////////////////////////////////////
/*                       PrincipalFusion                                */
//////////////////////////////////////////////////////////////////////////
// ����ת��  
void transMatrix( double *matrix, int m, int n );  
//�������
void multimatrix(double *a, double *b, int m, int n, int k, double *c);
// ��������  
void inverseMatrix( double *matrix, int n ); 

// ��ʵ�Գƾ��������ֵ�������������Ÿ�ȷ�  
bool eejcb( double a[], int n, double v[], double eps, int jt ); 

// ��������  
void linearStretch( unsigned char** pResult, int width, int height, int bandCount ); 
void linearStretch( double** pResult, int width, int height, int bandCount );
// ͳ�Ʋ��ξ�ֵ�����ͼ��Э�������
double* CovMatrix(const char* filepath);
  
// ������ֵ��С������������  
void sortEigenVector( double* eigenVector, double* covAfterEejcb,int nRastercount);  

// PCA�任  
//unsigned char** PCATransform( unsigned char **imgMatrix, double* eigenVector,int width,int height,int nRastercount);  

// PCA��任  
//unsigned char** inversePCA( unsigned char **imgMatrix, double* eigenVector,int width,int height,int n );  
// PCA�ں�  
//unsigned char** principalFusion(const char * filepath1,const char * filepath2) ;

//�����������ֵ����Сֵ
void max_min_value(double* array,int n, double& a,double& b);
#endif
















