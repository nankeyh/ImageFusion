/*
*程序功能介绍：实现图像融合处理操作，以及融合影像质量的客观评价。
*数据介绍：需要配准后（尺寸一致）不同分辨率（全色和多光谱）的8bit影像数据（程序提供GDAL库读写影像，支持库内提供的相关图像格式，本程序测试数据为tif影像）
*开发时间：2019.4.20
*开发者：YangHui
*可参考我的相关博客进行理解：https://blog.csdn.net/nanke_yh
*注意事项：本程序是全C/C++编写；
*          本程序重点测试的是相关的融合算法，影像质量评价相关函数在ImgQuality.cpp中，如若需要可自己调用；
*          本程序是根据8bit影像数据开发的，如果想改动，可根据自己需求将unsigned char数据类型进行相应的更改。
*本程序也是作者在参考许多资料的情况下基于自己的理解编写的，其中如若出现理解偏差或者错误，请查找更权威资料进行改正。也可联系作者进行交流qq：715643580
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

	//计时开始
	clock_t start = clock();
	//CorrelationFusion(MS,PAN,outfile[0]);

	principalFusion(MS,PAN,outfile[1]);

	//BroveyFusion(MS,PAN,outfile[2]);

	//IHSFusion(MS,PAN,outfile[3],0);

	//MultiplyFusion(MS,PAN,outfile[8],0);

	//HPFFusion(MS,PAN,outfile[7]);

	clock_t finish = clock(); //计时结束
	cout<<"总运行时间为："<<float(finish - start)/CLOCKS_PER_SEC<<"秒！"<<endl;
	system("pause");
	return 0;

}