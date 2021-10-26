#ifndef  _FUNCTION_H_
#define  _FUNCTION_H_

//#include "Macro.h"               
#include <iostream>              //载入头文件
#include <vector>                //矢量函数
#include <fstream>               //输出函数
#include <Windows.h>             //与计算机兼容函数
#include <string>                //字符串函数

//读取大影像
#include "../GDAL/gdal200/include/gdal.h"
#include "../GDAL/gdal200/include/gdal_priv.h"

using namespace std;

#define ColorLen    256         //影像灰度值范围2^8


#ifdef _WIN64                   //调用GDAL库
#pragma comment(lib,"../GDAL/gdal200/x64/lib/gdal_i.lib")
#else
#pragma comment(lib,"../GDAL/gdal200/x86/lib/gdal_i.lib")
#endif

unsigned char* OpenTif(const char *FilePath);

unsigned char** GdalOpenTif(const char *FilePath);

unsigned char** GDALreadImage(int &width,int &height,int &nRastercount,const char *filepath);


void WriteTif( unsigned char* pImageData,char *savepath,int width,int height);

void SaveImage(unsigned char** pImageData,char *savepath,int width,int height,int nRastercount);

bool OpenRaw(string FileName,unsigned char *ptImage,int width,int height);

bool SaveImgRaw(string FileName,unsigned char *ptNewImage,int width,int height);
bool SaveImgRaw(string FileName,double *ptNewImage,int width,int height);
void writeOutImgHdr(char* pDstImgFileName,int width,int height,int nChannels);

bool Depth_16to8(unsigned char *ptImage,int width,int height,string FileName);

char* findImageTypeGDAL(char *pImgFileName);

bool writeImageGDAL(char *pDstImgFileName,unsigned char *pImageData,int width,int height,int nChannels);

bool WriteTiff(const wchar_t * FileName,void *ptImage,size_t Width, size_t Height,  size_t Byte);
#endif