#include "Function.h"
#include <vector>
#include <atlstr.h>
#include "windows.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <time.h>
#include <string>
#include <assert.h>
#include <io.h>


using namespace std;
///////////////////////////////////////////////////////////////////////////////
/*                        图像的读取与保存                                   */
///////////////////////////////////////////////////////////////////////////////

unsigned char* OpenTif(const char *FilePath)
{
	GDALAllRegister();//注册、读取图像
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");//使之支持中文路径
	GDALDataset *poDataset = NULL;
	poDataset = (GDALDataset*)GDALOpen(FilePath,GA_ReadOnly);
	if(poDataset == NULL)
	{
		cout<<"无法打开影像！"<<endl;
		GDALDestroyDriverManager();
	}
	//获取图像数据的参数
	int width = poDataset->GetRasterXSize();
	int height = poDataset->GetRasterYSize();
	int nRastercount = poDataset->GetRasterCount();
	//开辟内存
	unsigned char *pImageData = new unsigned char[width * height];
	int bandList = {1};
	poDataset->RasterIO(GF_Read,0,0,width,height,pImageData,width,height,GDT_Byte,1,&bandList,0,0,0);
	//cout<<"单波段影像读入完成！"<<endl;
	//关闭GDAL库相关驱动和释放内存
	GDALClose(poDataset);
	return pImageData;
}

unsigned char** GdalOpenTif(const char *FilePath)
{
	GDALAllRegister();//注册、读取图像
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");//使之支持中文路径
	GDALDataset *poDataset = NULL;
	poDataset = (GDALDataset*)GDALOpen(FilePath,GA_ReadOnly);
	if(poDataset == NULL)
	{
		cout<<"无法打开影像！"<<endl;
		GDALDestroyDriverManager();
	}
	//获取图像数据的参数
	int width = poDataset->GetRasterXSize();
	int height = poDataset->GetRasterYSize();
	int nRastercount = poDataset->GetRasterCount();
	//开辟内存
	unsigned char **pImageData = new unsigned char *[nRastercount];
	if (nRastercount == 1)//单波段影像
	{
		int bandList = {1};
		*pImageData= new unsigned char[width*height];
		poDataset->RasterIO(GF_Read,0,0,width,height,*pImageData,width,height,GDT_Byte,1,&bandList,0,0,0);
		cout<<"单波段影像读入完成！"<<endl;
		//GDALClose(poDataset);
		//return pImageData;
	} 
	else if(nRastercount >= 3 )//多波段影像
	{
		for (int j=0;j< nRastercount;j++)
		{
			pImageData[j] =  new unsigned char[width*height];
		}
		for (int i = 1;i <= nRastercount;i++)
		{
			int bandList = {i};
			poDataset->RasterIO(GF_Read,0,0,width,height,pImageData[i-1],width,height,GDT_Byte,1,&bandList,0,0,0);
			//GDALRasterBand *pBand;
			//pBand = poDataset->GetRasterBand(i);
			//CPLErr error;
			//error = pBand->RasterIO(GF_Read,0,0,width,height,pImageData[i-1],width,height,GDT_Byte,0,0);
			//   if (error == CE_Failure)
			//   {
			//	cout<<"读取图像数据时失败！"<<endl;
			//	GDALDestroyDriverManager();
			//   }
		}
		cout<<"多光谱影像读入完成！"<<endl;	
	}
	//关闭GDAL库相关驱动和释放内存
	GDALClose(poDataset);
	return pImageData;
}

unsigned char** GDALreadImage(int &width,int &height,int &nRastercount,const char *filepath)
{
	GDALAllRegister();//注册、读取图像
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");//使之支持中文路径
	GDALDataset *poDataset = NULL;
	poDataset = (GDALDataset*)GDALOpen(filepath,GA_ReadOnly);
	if(poDataset == NULL)
	{
		cout<<"无法打开图像！"<<endl;
		GDALDestroyDriverManager();
	}
	//获取图像数据的参数
	width = poDataset->GetRasterXSize();
	height = poDataset->GetRasterYSize();
	nRastercount = poDataset->GetRasterCount();

	int i = 0,j=0;
	//开辟内存
	unsigned char **pImageData = new unsigned char *[nRastercount];
	if (nRastercount == 1)//单波段影像
	{
		int bandList = {1};
		*pImageData= new unsigned char[width*height];
		poDataset->RasterIO(GF_Read,0,0,width,height,*pImageData,width,height,GDT_Byte,1,&bandList,0,0,0);
		cout<<"全色影像读入完成！"<<endl;
		//GDALClose(poDataset);
		//return pImageData;
	} 
	else if(nRastercount >= 3 )//多波段影像
	{
		for (j=0;j< nRastercount;j++)
		{
			pImageData[j] =  new unsigned char[width*height];
		}
		for (i=1;i<= nRastercount;i++)
		{
			int bandList = {i};
			//int bandList = {nRastercount-i+1};//影像读入时，读入顺序反了，就需要在输出影像数据是颠倒波段顺序
			poDataset->RasterIO(GF_Read,0,0,width,height,pImageData[i-1],width,height,GDT_Byte,1,&bandList,0,0,0);
			//GDALRasterBand *pBand;
			//pBand = poDataset->GetRasterBand(i);
			//CPLErr error;
			//error = pBand->RasterIO(GF_Read,0,0,width,height,pImageData[i-1],width,height,GDT_Byte,0,0);
			//   if (error == CE_Failure)
			//   {
			//	cout<<"读取图像数据时失败！"<<endl;
			//	GDALDestroyDriverManager();
			//   }
		}
		cout<<"多光谱影像读入完成！"<<endl;	
	}
	//关闭GDAL库相关驱动和释放内存
	GDALClose(poDataset);poDataset = NULL;
	return pImageData;
}




void WriteTif( unsigned char* pImageData,char *savepath,int width,int height)  
{  
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");//使之支持中文路径
	char *GType = NULL;
	GType = findImageTypeGDAL(savepath);
	// 创建结果存储图像  
	if ( GType == NULL )  
	{  
		cout << "没有指定存储图像的格式，默认为GTiff" << endl;  
	}  
	const char *format ="GTiff";
	GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName( format ); 
	GDALDataset* WriteDataSet = poDriver->Create( savepath, width, height, 1, GDT_Byte, NULL );
	GDALRasterBand* pand = WriteDataSet->GetRasterBand(1);  
	pand->RasterIO( GF_Write, 0, 0, width, height, pImageData, width, height, GDT_Byte, 0, 0 ); 
	GDALClose( WriteDataSet );  
	cout<<"单波段影像保存完成！"<<endl;

}

void SaveImage( unsigned char** pImageData,char *savepath,int width,int height,int nRastercount)  
{  
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");//使之支持中文路径
	char *GType = NULL;
	GType = findImageTypeGDAL(savepath);
	// 创建结果存储图像  
	if ( GType == NULL )  
	{  
		cout << "没有指定存储图像的格式，默认为GTiff" << endl;  
	}  
	const char *format ="GTiff";
	GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName( format ); 
	if (nRastercount == 1)
	{
		GDALDataset* WriteDataSet = poDriver->Create( savepath, width, height, nRastercount, GDT_Byte, NULL );
		GDALRasterBand* pand = WriteDataSet->GetRasterBand(nRastercount);  
		pand->RasterIO( GF_Write, 0, 0, width, height, *pImageData, width, height, GDT_Byte, 0, 0 ); 
		GDALClose( WriteDataSet );  
		cout<<"单波段影像保存完成！"<<endl;
	} 
	else if(nRastercount >= 3)//控制输出影像波段数   
	{
		char **papszOptions = NULL;
		papszOptions = CSLSetNameValue( papszOptions,"PHOTOMETRIC", "RGB");//RGB模式,只是PS中初始打开时显示模式（若不设置用ps打开不显示成彩色图像）
		GDALDataset* WriteDataSet = poDriver->Create( savepath, width, height, nRastercount, GDT_Byte, papszOptions );//tif格式支持直接create，某些格式的图像不支持直接create   
		for ( int i = 0; i < nRastercount; i++ )  
		{  
			GDALRasterBand* pand = WriteDataSet->GetRasterBand( i + 1 );//这里写入影像波段的顺序是：外部内存中存储影像波段一一的写入。同时，这里RGB三波段的顺序是按1R2G3B来的。那么如果要对应，则在外部写入时注意！
			pand->RasterIO( GF_Write, 0, 0, width, height, pImageData[i], width, height, GDT_Byte, 0, 0 );  
		} 
		GDALClose( WriteDataSet );  
		cout<<"彩色影像保存完成！"<<endl;
	}
}  

bool OpenRaw(string FileName,unsigned char *ptImage,int width,int height)
{
	FILE *fp = NULL;
	if (0 != fopen_s(&fp,FileName.c_str(),"rb"))
	{
		MessageBox(NULL,_T("影像打开失败，请检查文件是否存在或者被破坏！"),_T("系统提示"),MB_OK|MB_ICONERROR);
		return false;
	}
	fread(ptImage,sizeof(unsigned char),width * height,fp);
	fclose(fp);
	return true;
}

bool SaveImgRaw(string FileName,unsigned char *ptNewImage,int width,int height)
{
	if (ptNewImage == NULL)
	{
		MessageBox(NULL,_T("写入内存失败！"),_T("系统提示"),MB_OK|MB_ICONERROR);
		return false;		
	}
	FILE *fp = NULL;
	fopen_s(&fp,FileName.c_str(),"wb");
	fwrite(ptNewImage,sizeof(unsigned char),width * height,fp);
	fclose(fp);
	cout<<"影像输出成功！"<<endl;
	return true;
}
bool SaveImgRaw(string FileName,double *ptNewImage,int width,int height)
{
	if (ptNewImage == NULL)
	{
		MessageBox(NULL,_T("写入内存失败！"),_T("系统提示"),MB_OK|MB_ICONERROR);
		return false;		
	}
	FILE *fp = NULL;
	fopen_s(&fp,FileName.c_str(),"wb");
	for (int i = 0;i < height;i++)
	{
		for (int j = 0; j < width;j++)
		{
			unsigned char tmp = (unsigned char)(ptNewImage[i*width+j]);
			fwrite(&tmp,sizeof(unsigned char),1,fp);
		}
	}
	fclose(fp);
	cout<<"影像输出成功！"<<endl;
	return true;
}
void writeOutImgHdr(char* pDstImgFileName,int width,int height,int nChannels)
{
	//生成头文件
	string path_out_raw;
	string path_out_hdr;
	char path_out_raw_c[_MAX_PATH]={NULL};
	strcpy(path_out_raw_c,pDstImgFileName);
	path_out_raw=path_out_raw_c;
	int pos=path_out_raw.find(".raw");
	if (pos!=-1)
	path_out_hdr=path_out_raw.replace(path_out_raw.find(".raw"),4,".hdr");
	else
	path_out_hdr=strcat(path_out_raw_c,".hdr");

	FILE * fp_hdr=fopen(path_out_hdr.c_str(),"w+t");
	fprintf(fp_hdr,"ENVI\ndescription = {\n  File Imported into ENVI.}\nsamples =  %d\nlines   =  %d\nbands   = %d\nheader offset = 0\nfile type = ENVI Standard\ndata type = 12\ninterleave = bsq\nsensor type = Unknown\nbyte order = 0\nwavelength units = Unknown",width,height,nChannels);
	fclose(fp_hdr);
}
bool Depth_16to8(unsigned char *ptImage,int width,int height,string FileName)
{
	unsigned char *img = new unsigned char[width * height];
	memset(img,0,sizeof(unsigned char) * width * height);
	for (int i = 0;i < height;i++)
	{
		for (int j = 0;j < width;j++)
		{
			img[i * width + j] = (unsigned char)(ptImage[i * width + j]/257); // 65535/255 = 257
		}
	}
	FILE *fp = NULL;
	fopen_s(&fp,FileName.c_str(),"wb");
	fwrite(img,sizeof(unsigned char),width*height,fp);
	fclose(fp);
	//cout<<"影像输出成功！"<<endl;
	delete []img;
	return true;
}

/******************************************************************************
函数名：
	findImageTypeGDAL
功能：
	根据文件后缀判断输入输出文件为GDAL图像格式还是裸数据格式
参数：
	char *pImgFileName       - 文件路径名称
说明：**************************************************************************/
char* findImageTypeGDAL(char *pImgFileName)
{
	char *Gtype = NULL;
	char fileExtension[_MAX_PATH]={NULL};
	strcpy(fileExtension,pImgFileName);
	if      ((strstr(fileExtension,".tif")!=NULL))  Gtype = "GTiff";
	else if ((strstr(fileExtension,".tiff")!=NULL)) Gtype = "GTiff";
	else if ((strstr(fileExtension,".jpg")!=NULL))  Gtype = "JPEG";
	else if ((strstr(fileExtension,".bmp")!=NULL))  Gtype = "BMP";
	else if ((strstr(fileExtension,".png")!=NULL))  Gtype = "PNG";  
	else if ((strstr(fileExtension,".gif")!=NULL))  Gtype = "GIF";
	else Gtype = NULL;
	//else Gtype = "GDALRaw";
	return Gtype;
}

/******************************************************************************
函数名：
	writeImageGDAL
功能：
	保存数据为GDAL图像格式(相当于将裸数据转换成tiff格式)
参数：
	unsigned char *pImageData        - 指向图像数据的指针
	int width,int height        - 图像宽度、高度、波段数
	int nChannels               - 图像通道，一般为1或3。1代表灰度图像，3代表RGB图像。
	char *pDstImgFileName       - 输出文件路径名称
说明：******************************************************************************/
bool writeImageGDAL(char *pDstImgFileName,unsigned char *pImageData,int width,int height,int nChannels)
{
	char headerfile[_MAX_PATH];	
	sprintf(headerfile,"tiffhead\\width%d_height%d_band%d.head",width,height,nChannels);
	if(access("tiffhead",0)==-1)
	{   //不存在就创建
		mkdir("tiffhead");//文件夹名
	}
	if(access(headerfile,0)!=-1)
	{   //判断tiff头文件是否存在
		long headersize = filelength(open(headerfile, 0x0100)); //open file for read,get length of file
		char*headbuffer=new char[headersize]();

		FILE *fp_head= fopen(headerfile,"rb+"); 
		fread(headbuffer,sizeof(char),headersize,fp_head);//读取头文件
		fclose(fp_head);

		//写raw的方式写tiff
		FILE *fp_out=fopen(pDstImgFileName,"wb+");
		if (fp_out==NULL)
			return false;
		fwrite(headbuffer,sizeof(char),headersize,fp_out);
		fclose(fp_out);
		fp_out= fopen(pDstImgFileName,"ab+");
		fwrite(pImageData,sizeof(unsigned char),height*width*nChannels,fp_out); 
		fclose(fp_out);
		
		delete []headbuffer;
	}
	else
	{
		assert ( !(pDstImgFileName == NULL || pImageData == NULL || width <1 || height < 1 || nChannels < 1));

		GDALAllRegister();
		CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");//使之支持中文路径
		char *GType = NULL;
		GType = findImageTypeGDAL(pDstImgFileName);
		if (GType == NULL)	
		{ 
			return false;
		}
		//const char *pszFormat=GType;
		const char *pszFormat ="GTiff";
		GDALDriver *poDriver=GetGDALDriverManager()->GetDriverByName(pszFormat);
		if( poDriver == NULL ) 
			return false; 

		char **papszOptions = NULL;
		if (nChannels==3)
		{
			papszOptions = CSLSetNameValue( papszOptions, "PHOTOMETRIC", "RGB" );//三通道图像设置为RGB模式，若不设置用ps打开不显示成彩色图像
		}
		GDALDataset *WriteDataSet = poDriver->Create(pDstImgFileName, width,height,nChannels,GDT_Byte,papszOptions);//tif格式支持直接create，某些格式的图像不支持直接create
		WriteDataSet->RasterIO(GF_Write,0,0,width,height,pImageData,width,height,GDT_Byte,nChannels,NULL,0,0,0);
		//保存图像的相关信息
		//GDALClose(WriteDataSet);
		delete WriteDataSet;
		WriteDataSet=NULL;

		//生成header文件,首先计算出头文件大小
		HANDLE handle = CreateFileA(pDstImgFileName, FILE_READ_EA, FILE_SHARE_READ, 0, OPEN_EXISTING, 0, 0);
		if (handle != INVALID_HANDLE_VALUE)
		{
			unsigned long tmp1=GetFileSize(handle, NULL);
			unsigned long tmp2 =width*height*nChannels*sizeof(char);//图像裸数据大小
			unsigned long headsize = tmp1-tmp2; //tif图像大小减去裸数据大小获得的就是图像头信息
			CloseHandle(handle);

			char*headbuffer=new char[headsize]();
			FILE *fp= fopen(pDstImgFileName,"rb+"); 
			fread(headbuffer,sizeof(char),headsize,fp);//读取图像头信息
			fclose(fp);

			FILE *fp_head=fopen(headerfile,"wb");
			fwrite(headbuffer,sizeof(char),headsize,fp_head); 
			fclose(fp_head);
			delete []headbuffer;
		}
	}
	cout<<"影像输出完成！"<<endl;
	return true;
}

/************************************************************************************************************
@func：对数据指针，保存成tiff文件格式 （仅能保存单波段图像，可实现图像波段的分离，多波段图像保存问题需再研究）
@param：void *ptImage，保存数据的影像指针
			size_t Height 影像的高度
			size_t Width  影像的宽度
			size_t Byte   影像的位数
@return： 如果正确的话，返回true，不然的返回false
图像文件头Image File Header(IFH),图像文件目录Image File Directory(IFD)和目录项Directory Entry(DE)。
说明：*******************************************************************************************************/
bool WriteTiff(const wchar_t * FileName,void *ptImage, size_t Height, size_t Width, size_t Byte)
{
	if (NULL==ptImage)
	{
		cout << "image open failed!" << endl;
		return false;
	}
	//CreateFile存在两种情况：CreateFileA和CreateFileW，在此用CreateFileW直接注明
	HANDLE handle = CreateFileW(FileName, GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	if (INVALID_HANDLE_VALUE==handle)
	{
		cout << "Create  FileName  Failed！" << endl; 
		return false;
	}
	//写IFH 8字节
	//4949
	unsigned char temp = 0x4949;
	unsigned char *Order = (unsigned char*)&temp;
	DWORD WriteSize = 0;
	WriteFile(handle, Order, 2, &WriteSize, NULL);
	if (2!=WriteSize)
	{
		cout << "Write the Order failed" << endl;
		return false;
	}
	//version
	temp = 0x002A;
	WriteFile(handle, Order, 2, &WriteSize, NULL);
	if (2!=WriteSize)
	{
		cout << "Write the version info failed" << endl;
		return false;
	}
	DWORD Length = 8+Byte*Height*Width; //第一个IFD的位置  --等于文件的指针8+ptimage的大小
	// Write the tiff file offset
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	if (4!=WriteSize)
	{
		cout << "Write the Offset Failed!" << endl;
		return false;
	}
	//写入影像的数据
	WriteFile(handle, ptImage, Byte*Width*Height, &WriteSize, NULL);
	//写入IFD文件
	//DE数量
	//宽，高，分辨率，像素位置，是否压缩，反色，自己的信息 一共是七个信息
	temp = 0x07;
	WriteFile(handle, &temp, 2, &WriteSize, NULL);
	if (2!=WriteSize)
	{
		cout << "Write the DE Count failed!" << endl;
		return false;
	}

	unsigned char tag = 0;
	unsigned char type = 0;
	//unsigned Length = 0;
	unsigned ValueOffset = 0;
	//写入影像的宽
	tag = 0x0100;
	type = 0x0003;
	Length = 0x01;
	ValueOffset = Width;
	WriteFile(handle, &tag, 2, &WriteSize, NULL);
	WriteFile(handle, &type, 2, &WriteSize, NULL);
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	WriteFile(handle, &ValueOffset, 4, &WriteSize, NULL);
	//写入影像的高
	tag = 0x0101;
	type = 0x0003;
	Length = 0x01;
	ValueOffset = Height;
	WriteFile(handle, &tag, 2, &WriteSize, NULL);
	WriteFile(handle, &type, 2, &WriteSize, NULL);
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	WriteFile(handle, &ValueOffset, 4, &WriteSize, NULL);
	//写入影像的分辨率
	tag = 0x0102;
	type = 0x0003;
	Length = 0x01;       //类型数据个数控制真彩色（>2即为彩色）
	ValueOffset = Byte*8;//位深（值=1为单色，=4为16色，=8为256色）
	WriteFile(handle, &tag, 2, &WriteSize, NULL);
	WriteFile(handle, &type, 2, &WriteSize, NULL);
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	WriteFile(handle, &ValueOffset, 4, &WriteSize, NULL);
	//写入是否压缩(这一项意义不大，没有也是以不压缩输出)
	tag = 0x0103;
	type = 0x0003;
	Length = 0x01;
	ValueOffset = 0x08;//（05表示压缩，否则不压缩）
	WriteFile(handle, &tag, 2, &WriteSize, NULL);
	WriteFile(handle, &type, 2, &WriteSize, NULL);
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	WriteFile(handle, &ValueOffset, 4, &WriteSize, NULL);
	//写入影像是反色的
	/*不加这一个信息的话，输出文件将是原始图像的反色情况*/
	tag = 0x0106;
	type = 0x0003;
	Length = 0x01;
	ValueOffset = 0x01;//（01表示反色，否则不反色）
	WriteFile(handle, &tag, 2, &WriteSize, NULL);
	WriteFile(handle, &type, 2, &WriteSize, NULL);
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	WriteFile(handle, &ValueOffset, 4, &WriteSize, NULL);
	//写入数据所在的位置
	tag = 0x0111;
	type = 0x0004;
	Length = 0x01;
	ValueOffset = 0x08;
	WriteFile(handle, &tag, 2, &WriteSize, NULL);
	WriteFile(handle, &type, 2, &WriteSize, NULL);
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	WriteFile(handle, &ValueOffset, 4, &WriteSize, NULL);
	//写入自己的信息
	tag = 0x013B;
	type = 0x0002;
	Length = 8;//信息字节数+1
	//8  前8个字节;宽*高*字节数 像素的大小;2  Directory Entry Count占用的2字节;DE的数量7*12  每个IFD所占的字节数 ;4  下一个tag的位置Offset to next IFD;
	ValueOffset = 8+Height*Width*Byte+2+12*7+4;
	WriteFile(handle, &tag, 2, &WriteSize, NULL);
	WriteFile(handle, &type, 2, &WriteSize, NULL);
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	WriteFile(handle, &ValueOffset, 4, &WriteSize, NULL);
	//写入文件尾
	Length = 0x00;
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	//写入自己的信息
	char *text="YANGHUI";
	WriteFile(handle,text,8,&WriteSize,NULL);
	cout<<"单波段影像输出完成！"<<endl;
	return true;
}

 



























