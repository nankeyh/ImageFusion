/*
*�����ܽ��ܣ�ʵ��ͼ���ںϴ���������Լ��ں�Ӱ�������Ŀ͹����ۡ�
*���ݽ��ܣ���Ҫ��׼�󣨳ߴ�һ�£���ͬ�ֱ��ʣ�ȫɫ�Ͷ���ף���8bitӰ�����ݣ������ṩGDAL���дӰ��֧�ֿ����ṩ�����ͼ���ʽ���������������ΪtifӰ��
*����ʱ�䣺2019.4.20
*�����ߣ�YangHui
*�ɲο��ҵ���ز��ͽ�����⣺https://blog.csdn.net/nanke_yh
*ע�������������ȫC/C++��д��
*          �������ص���Ե�����ص��ں��㷨��Ӱ������������غ�����ImgQuality.cpp�У�������Ҫ���Լ����ã�
*          �������Ǹ���8bitӰ�����ݿ����ģ������Ķ����ɸ����Լ�����unsigned char�������ͽ�����Ӧ�ĸ��ġ�
*������Ҳ�������ڲο�������ϵ�����»����Լ�������д�ģ����������������ƫ����ߴ�������Ҹ�Ȩ�����Ͻ��и�����Ҳ����ϵ���߽��н���qq��715643580
*/
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <Windows.h>
#include <ctime>
//#include "Macro.h"
#include "Function.h"
#include "Fusion.h"
#include <stdlib.h>
#include "ImgQuality.h"
using namespace std;

int main(int argc, char* argv[])
{
	char *PAN = "data/spot.tif";
	char *MS = "data/tm.tif";

	char *outfile[10] = {"result/ImgFusion_correlation.tif","result/ImgFusion_PCA.tif","result/ImgFusion_Brovey.tif",
	                   "result/ImgFusion_IHS0.tif","result/ImgFusion_IHS1.tif","result/ImgFusion_IHS2.tif","result/ImgFusion_IHS3.tif",
					   "result/ImgFusion_HPF.tif","result/ImgFusion_Multiply0.tif","result/ImgFusion_Multiply1.tif"}; 

	//��ʱ��ʼ
	clock_t start = clock();
	//CorrelationFusion(MS,PAN,outfile[0]);

	principalFusion(MS,PAN,outfile[1]);

	//BroveyFusion(MS,PAN,outfile[2]);

	//IHSFusion(MS,PAN,outfile[3],0);

	//MultiplyFusion(MS,PAN,outfile[8],0);

	//HPFFusion(MS,PAN,outfile[7]);

	clock_t finish = clock(); //��ʱ����
	cout<<"������ʱ��Ϊ��"<<float(finish - start)/CLOCKS_PER_SEC<<"�룡"<<endl;
	system("pause");
	return 0;

}