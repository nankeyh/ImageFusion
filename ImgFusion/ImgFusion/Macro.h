#pragma once                                  //��֤�궨������ֻ��һ��
#ifndef  __IFMACRO_H__			              
#define  __IFMACRO_H__			              

#include <cmath>
//�������ã�unsigned char	                  

/*�곣��*/						               
#define ImageBits   8                        //Ӱ��λ��

#define ColorLen    256                     //Ӱ��Ҷ�ֵ��Χ2^8
//const int ImgHeight = 331;                  //��׼��Ӱ��ĸ�
//const int ImgWidth = 331;                   //��׼��Ӱ��Ŀ�
//const int ImgHeight = 600;                  //��׼��Ӱ��ĸ�
//const int ImgWidth = 400;                   //��׼��Ӱ��Ŀ�

//const int  Num = ImgWidth*ImgHeight;        //��ߺ궨����Ҫ���ڼ���width*height������ֻ����width*height

const int FilterWidth = 3;		              //�˲����ڵĳߴ�

const double PI = acos(-1.0);                 //Բ���ʦ�ֵ3.1415926


#endif
