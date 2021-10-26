#ifndef  _IMGQUALITY_H_
#define  _IMGQUALITY_H_

#include "Function.h"
#include "Macro.h"           
#include <iostream>          
#include <vector>            
#include <fstream>           
#include <Windows.h>         
#include <string>            

using namespace std;

/***********************计算 calculation*********************************/
/************************************************************************/
/*                      1、均值计算（单片23434*23619）                  */
/************************************************************************/
//bool CalAverage(unsigned char *ptImage,double &ImgAver); //取址符直接获得均值
double CalAverage(unsigned char *ptImage,int width,int height); //返回均值

/************************************************************************/
/*                      2、标准差计算                                   */
/************************************************************************/
double CalStandardDevi(unsigned char *ptImage,int width,int height);

/************************************************************************/
/*                      3、空间频率计算                                 */
/************************************************************************/
double CalSpatialFreq(unsigned char *ptImage,int width,int height);

/************************************************************************/
/*                      4、平均梯度计算                                 */
/************************************************************************/
double CalAverGrad(unsigned char *ptImage,int width,int height);

/************************************************************************/
/*                      5、信息熵计算                                   */
/************************************************************************/
double CalEntropy(unsigned char *ptImage,int width,int height);

/************************************************************************/
/*                      6、相关系数计算                                 */
/************************************************************************/
double CalCorrelationCoef(unsigned char *ptImage,unsigned char *refImage,int width,int height);

/************************************************************************/
/*                      7、偏差指数计算                                 */
/************************************************************************/
double CalDeviationIndex(unsigned char *dstImage,unsigned char *refImage,int width,int height);

/************************************************************************/
/*                      8、扭曲度计算                                   */
/************************************************************************/
double CalDistortion(unsigned char *dstImage,unsigned char *refImage,int width,int height);










#endif





















