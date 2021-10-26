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
// 矩阵转置  
void transMatrix( double *matrix, int m, int n );  
//矩阵相乘
void multimatrix(double *a, double *b, int m, int n, int k, double *c);
// 矩阵求逆  
void inverseMatrix( double *matrix, int n ); 

// 求实对称矩阵的特征值及特征向量的雅格比法  
bool eejcb( double a[], int n, double v[], double eps, int jt ); 

// 线性拉伸  
void linearStretch( unsigned char** pResult, int width, int height, int bandCount ); 
void linearStretch( double** pResult, int width, int height, int bandCount );
// 统计波段均值并求得图像协方差矩阵
double* CovMatrix(const char* filepath);
  
// 按特征值大小排序特征向量  
void sortEigenVector( double* eigenVector, double* covAfterEejcb,int nRastercount);  

// PCA变换  
//unsigned char** PCATransform( unsigned char **imgMatrix, double* eigenVector,int width,int height,int nRastercount);  

// PCA逆变换  
//unsigned char** inversePCA( unsigned char **imgMatrix, double* eigenVector,int width,int height,int n );  
// PCA融合  
//unsigned char** principalFusion(const char * filepath1,const char * filepath2) ;

//求多个数的最大值和最小值
void max_min_value(double* array,int n, double& a,double& b);
#endif
















