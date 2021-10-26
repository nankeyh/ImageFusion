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

/***********************���� calculation*********************************/
/************************************************************************/
/*                      1����ֵ���㣨��Ƭ23434*23619��                  */
/************************************************************************/
//bool CalAverage(unsigned char *ptImage,double &ImgAver); //ȡַ��ֱ�ӻ�þ�ֵ
double CalAverage(unsigned char *ptImage,int width,int height); //���ؾ�ֵ

/************************************************************************/
/*                      2����׼�����                                   */
/************************************************************************/
double CalStandardDevi(unsigned char *ptImage,int width,int height);

/************************************************************************/
/*                      3���ռ�Ƶ�ʼ���                                 */
/************************************************************************/
double CalSpatialFreq(unsigned char *ptImage,int width,int height);

/************************************************************************/
/*                      4��ƽ���ݶȼ���                                 */
/************************************************************************/
double CalAverGrad(unsigned char *ptImage,int width,int height);

/************************************************************************/
/*                      5����Ϣ�ؼ���                                   */
/************************************************************************/
double CalEntropy(unsigned char *ptImage,int width,int height);

/************************************************************************/
/*                      6�����ϵ������                                 */
/************************************************************************/
double CalCorrelationCoef(unsigned char *ptImage,unsigned char *refImage,int width,int height);

/************************************************************************/
/*                      7��ƫ��ָ������                                 */
/************************************************************************/
double CalDeviationIndex(unsigned char *dstImage,unsigned char *refImage,int width,int height);

/************************************************************************/
/*                      8��Ť���ȼ���                                   */
/************************************************************************/
double CalDistortion(unsigned char *dstImage,unsigned char *refImage,int width,int height);










#endif





















