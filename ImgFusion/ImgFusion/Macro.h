#pragma once                                  //保证宏定义有且只有一个
#ifndef  __IFMACRO_H__			              
#define  __IFMACRO_H__			              

#include <cmath>
//数据类用：unsigned char	                  

/*宏常量*/						               
#define ImageBits   8                        //影像位深

#define ColorLen    256                     //影像灰度值范围2^8
//const int ImgHeight = 331;                  //配准后影像的高
//const int ImgWidth = 331;                   //配准后影像的宽
//const int ImgHeight = 600;                  //配准后影像的高
//const int ImgWidth = 400;                   //配准后影像的宽

//const int  Num = ImgWidth*ImgHeight;        //宽高宏定义主要用于计算width*height，后续只调用width*height

const int FilterWidth = 3;		              //滤波窗口的尺寸

const double PI = acos(-1.0);                 //圆周率π值3.1415926


#endif
