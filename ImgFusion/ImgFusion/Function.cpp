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
/*                        ͼ��Ķ�ȡ�뱣��                                   */
///////////////////////////////////////////////////////////////////////////////

unsigned char* OpenTif(const char *FilePath)
{
	GDALAllRegister();//ע�ᡢ��ȡͼ��
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");//ʹ֧֮������·��
	GDALDataset *poDataset = NULL;
	poDataset = (GDALDataset*)GDALOpen(FilePath,GA_ReadOnly);
	if(poDataset == NULL)
	{
		cout<<"�޷���Ӱ��"<<endl;
		GDALDestroyDriverManager();
	}
	//��ȡͼ�����ݵĲ���
	int width = poDataset->GetRasterXSize();
	int height = poDataset->GetRasterYSize();
	int nRastercount = poDataset->GetRasterCount();
	//�����ڴ�
	unsigned char *pImageData = new unsigned char[width * height];
	int bandList = {1};
	poDataset->RasterIO(GF_Read,0,0,width,height,pImageData,width,height,GDT_Byte,1,&bandList,0,0,0);
	//cout<<"������Ӱ�������ɣ�"<<endl;
	//�ر�GDAL������������ͷ��ڴ�
	GDALClose(poDataset);
	return pImageData;
}

unsigned char** GdalOpenTif(const char *FilePath)
{
	GDALAllRegister();//ע�ᡢ��ȡͼ��
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");//ʹ֧֮������·��
	GDALDataset *poDataset = NULL;
	poDataset = (GDALDataset*)GDALOpen(FilePath,GA_ReadOnly);
	if(poDataset == NULL)
	{
		cout<<"�޷���Ӱ��"<<endl;
		GDALDestroyDriverManager();
	}
	//��ȡͼ�����ݵĲ���
	int width = poDataset->GetRasterXSize();
	int height = poDataset->GetRasterYSize();
	int nRastercount = poDataset->GetRasterCount();
	//�����ڴ�
	unsigned char **pImageData = new unsigned char *[nRastercount];
	if (nRastercount == 1)//������Ӱ��
	{
		int bandList = {1};
		*pImageData= new unsigned char[width*height];
		poDataset->RasterIO(GF_Read,0,0,width,height,*pImageData,width,height,GDT_Byte,1,&bandList,0,0,0);
		cout<<"������Ӱ�������ɣ�"<<endl;
		//GDALClose(poDataset);
		//return pImageData;
	} 
	else if(nRastercount >= 3 )//�ನ��Ӱ��
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
			//	cout<<"��ȡͼ������ʱʧ�ܣ�"<<endl;
			//	GDALDestroyDriverManager();
			//   }
		}
		cout<<"�����Ӱ�������ɣ�"<<endl;	
	}
	//�ر�GDAL������������ͷ��ڴ�
	GDALClose(poDataset);
	return pImageData;
}

unsigned char** GDALreadImage(int &width,int &height,int &nRastercount,const char *filepath)
{
	GDALAllRegister();//ע�ᡢ��ȡͼ��
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");//ʹ֧֮������·��
	GDALDataset *poDataset = NULL;
	poDataset = (GDALDataset*)GDALOpen(filepath,GA_ReadOnly);
	if(poDataset == NULL)
	{
		cout<<"�޷���ͼ��"<<endl;
		GDALDestroyDriverManager();
	}
	//��ȡͼ�����ݵĲ���
	width = poDataset->GetRasterXSize();
	height = poDataset->GetRasterYSize();
	nRastercount = poDataset->GetRasterCount();

	int i = 0,j=0;
	//�����ڴ�
	unsigned char **pImageData = new unsigned char *[nRastercount];
	if (nRastercount == 1)//������Ӱ��
	{
		int bandList = {1};
		*pImageData= new unsigned char[width*height];
		poDataset->RasterIO(GF_Read,0,0,width,height,*pImageData,width,height,GDT_Byte,1,&bandList,0,0,0);
		cout<<"ȫɫӰ�������ɣ�"<<endl;
		//GDALClose(poDataset);
		//return pImageData;
	} 
	else if(nRastercount >= 3 )//�ನ��Ӱ��
	{
		for (j=0;j< nRastercount;j++)
		{
			pImageData[j] =  new unsigned char[width*height];
		}
		for (i=1;i<= nRastercount;i++)
		{
			int bandList = {i};
			//int bandList = {nRastercount-i+1};//Ӱ�����ʱ������˳���ˣ�����Ҫ�����Ӱ�������ǵߵ�����˳��
			poDataset->RasterIO(GF_Read,0,0,width,height,pImageData[i-1],width,height,GDT_Byte,1,&bandList,0,0,0);
			//GDALRasterBand *pBand;
			//pBand = poDataset->GetRasterBand(i);
			//CPLErr error;
			//error = pBand->RasterIO(GF_Read,0,0,width,height,pImageData[i-1],width,height,GDT_Byte,0,0);
			//   if (error == CE_Failure)
			//   {
			//	cout<<"��ȡͼ������ʱʧ�ܣ�"<<endl;
			//	GDALDestroyDriverManager();
			//   }
		}
		cout<<"�����Ӱ�������ɣ�"<<endl;	
	}
	//�ر�GDAL������������ͷ��ڴ�
	GDALClose(poDataset);poDataset = NULL;
	return pImageData;
}




void WriteTif( unsigned char* pImageData,char *savepath,int width,int height)  
{  
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");//ʹ֧֮������·��
	char *GType = NULL;
	GType = findImageTypeGDAL(savepath);
	// ��������洢ͼ��  
	if ( GType == NULL )  
	{  
		cout << "û��ָ���洢ͼ��ĸ�ʽ��Ĭ��ΪGTiff" << endl;  
	}  
	const char *format ="GTiff";
	GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName( format ); 
	GDALDataset* WriteDataSet = poDriver->Create( savepath, width, height, 1, GDT_Byte, NULL );
	GDALRasterBand* pand = WriteDataSet->GetRasterBand(1);  
	pand->RasterIO( GF_Write, 0, 0, width, height, pImageData, width, height, GDT_Byte, 0, 0 ); 
	GDALClose( WriteDataSet );  
	cout<<"������Ӱ�񱣴���ɣ�"<<endl;

}

void SaveImage( unsigned char** pImageData,char *savepath,int width,int height,int nRastercount)  
{  
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");//ʹ֧֮������·��
	char *GType = NULL;
	GType = findImageTypeGDAL(savepath);
	// ��������洢ͼ��  
	if ( GType == NULL )  
	{  
		cout << "û��ָ���洢ͼ��ĸ�ʽ��Ĭ��ΪGTiff" << endl;  
	}  
	const char *format ="GTiff";
	GDALDriver *poDriver = GetGDALDriverManager()->GetDriverByName( format ); 
	if (nRastercount == 1)
	{
		GDALDataset* WriteDataSet = poDriver->Create( savepath, width, height, nRastercount, GDT_Byte, NULL );
		GDALRasterBand* pand = WriteDataSet->GetRasterBand(nRastercount);  
		pand->RasterIO( GF_Write, 0, 0, width, height, *pImageData, width, height, GDT_Byte, 0, 0 ); 
		GDALClose( WriteDataSet );  
		cout<<"������Ӱ�񱣴���ɣ�"<<endl;
	} 
	else if(nRastercount >= 3)//�������Ӱ�񲨶���   
	{
		char **papszOptions = NULL;
		papszOptions = CSLSetNameValue( papszOptions,"PHOTOMETRIC", "RGB");//RGBģʽ,ֻ��PS�г�ʼ��ʱ��ʾģʽ������������ps�򿪲���ʾ�ɲ�ɫͼ��
		GDALDataset* WriteDataSet = poDriver->Create( savepath, width, height, nRastercount, GDT_Byte, papszOptions );//tif��ʽ֧��ֱ��create��ĳЩ��ʽ��ͼ��֧��ֱ��create   
		for ( int i = 0; i < nRastercount; i++ )  
		{  
			GDALRasterBand* pand = WriteDataSet->GetRasterBand( i + 1 );//����д��Ӱ�񲨶ε�˳���ǣ��ⲿ�ڴ��д洢Ӱ�񲨶�һһ��д�롣ͬʱ������RGB�����ε�˳���ǰ�1R2G3B���ġ���ô���Ҫ��Ӧ�������ⲿд��ʱע�⣡
			pand->RasterIO( GF_Write, 0, 0, width, height, pImageData[i], width, height, GDT_Byte, 0, 0 );  
		} 
		GDALClose( WriteDataSet );  
		cout<<"��ɫӰ�񱣴���ɣ�"<<endl;
	}
}  

bool OpenRaw(string FileName,unsigned char *ptImage,int width,int height)
{
	FILE *fp = NULL;
	if (0 != fopen_s(&fp,FileName.c_str(),"rb"))
	{
		MessageBox(NULL,_T("Ӱ���ʧ�ܣ������ļ��Ƿ���ڻ��߱��ƻ���"),_T("ϵͳ��ʾ"),MB_OK|MB_ICONERROR);
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
		MessageBox(NULL,_T("д���ڴ�ʧ�ܣ�"),_T("ϵͳ��ʾ"),MB_OK|MB_ICONERROR);
		return false;		
	}
	FILE *fp = NULL;
	fopen_s(&fp,FileName.c_str(),"wb");
	fwrite(ptNewImage,sizeof(unsigned char),width * height,fp);
	fclose(fp);
	cout<<"Ӱ������ɹ���"<<endl;
	return true;
}
bool SaveImgRaw(string FileName,double *ptNewImage,int width,int height)
{
	if (ptNewImage == NULL)
	{
		MessageBox(NULL,_T("д���ڴ�ʧ�ܣ�"),_T("ϵͳ��ʾ"),MB_OK|MB_ICONERROR);
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
	cout<<"Ӱ������ɹ���"<<endl;
	return true;
}
void writeOutImgHdr(char* pDstImgFileName,int width,int height,int nChannels)
{
	//����ͷ�ļ�
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
	//cout<<"Ӱ������ɹ���"<<endl;
	delete []img;
	return true;
}

/******************************************************************************
��������
	findImageTypeGDAL
���ܣ�
	�����ļ���׺�ж���������ļ�ΪGDALͼ���ʽ���������ݸ�ʽ
������
	char *pImgFileName       - �ļ�·������
˵����**************************************************************************/
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
��������
	writeImageGDAL
���ܣ�
	��������ΪGDALͼ���ʽ(�൱�ڽ�������ת����tiff��ʽ)
������
	unsigned char *pImageData        - ָ��ͼ�����ݵ�ָ��
	int width,int height        - ͼ���ȡ��߶ȡ�������
	int nChannels               - ͼ��ͨ����һ��Ϊ1��3��1����Ҷ�ͼ��3����RGBͼ��
	char *pDstImgFileName       - ����ļ�·������
˵����******************************************************************************/
bool writeImageGDAL(char *pDstImgFileName,unsigned char *pImageData,int width,int height,int nChannels)
{
	char headerfile[_MAX_PATH];	
	sprintf(headerfile,"tiffhead\\width%d_height%d_band%d.head",width,height,nChannels);
	if(access("tiffhead",0)==-1)
	{   //�����ھʹ���
		mkdir("tiffhead");//�ļ�����
	}
	if(access(headerfile,0)!=-1)
	{   //�ж�tiffͷ�ļ��Ƿ����
		long headersize = filelength(open(headerfile, 0x0100)); //open file for read,get length of file
		char*headbuffer=new char[headersize]();

		FILE *fp_head= fopen(headerfile,"rb+"); 
		fread(headbuffer,sizeof(char),headersize,fp_head);//��ȡͷ�ļ�
		fclose(fp_head);

		//дraw�ķ�ʽдtiff
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
		CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");//ʹ֧֮������·��
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
			papszOptions = CSLSetNameValue( papszOptions, "PHOTOMETRIC", "RGB" );//��ͨ��ͼ������ΪRGBģʽ������������ps�򿪲���ʾ�ɲ�ɫͼ��
		}
		GDALDataset *WriteDataSet = poDriver->Create(pDstImgFileName, width,height,nChannels,GDT_Byte,papszOptions);//tif��ʽ֧��ֱ��create��ĳЩ��ʽ��ͼ��֧��ֱ��create
		WriteDataSet->RasterIO(GF_Write,0,0,width,height,pImageData,width,height,GDT_Byte,nChannels,NULL,0,0,0);
		//����ͼ��������Ϣ
		//GDALClose(WriteDataSet);
		delete WriteDataSet;
		WriteDataSet=NULL;

		//����header�ļ�,���ȼ����ͷ�ļ���С
		HANDLE handle = CreateFileA(pDstImgFileName, FILE_READ_EA, FILE_SHARE_READ, 0, OPEN_EXISTING, 0, 0);
		if (handle != INVALID_HANDLE_VALUE)
		{
			unsigned long tmp1=GetFileSize(handle, NULL);
			unsigned long tmp2 =width*height*nChannels*sizeof(char);//ͼ�������ݴ�С
			unsigned long headsize = tmp1-tmp2; //tifͼ���С��ȥ�����ݴ�С��õľ���ͼ��ͷ��Ϣ
			CloseHandle(handle);

			char*headbuffer=new char[headsize]();
			FILE *fp= fopen(pDstImgFileName,"rb+"); 
			fread(headbuffer,sizeof(char),headsize,fp);//��ȡͼ��ͷ��Ϣ
			fclose(fp);

			FILE *fp_head=fopen(headerfile,"wb");
			fwrite(headbuffer,sizeof(char),headsize,fp_head); 
			fclose(fp_head);
			delete []headbuffer;
		}
	}
	cout<<"Ӱ�������ɣ�"<<endl;
	return true;
}

/************************************************************************************************************
@func��������ָ�룬�����tiff�ļ���ʽ �����ܱ��浥����ͼ�񣬿�ʵ��ͼ�񲨶εķ��룬�ನ��ͼ�񱣴����������о���
@param��void *ptImage���������ݵ�Ӱ��ָ��
			size_t Height Ӱ��ĸ߶�
			size_t Width  Ӱ��Ŀ��
			size_t Byte   Ӱ���λ��
@return�� �����ȷ�Ļ�������true����Ȼ�ķ���false
ͼ���ļ�ͷImage File Header(IFH),ͼ���ļ�Ŀ¼Image File Directory(IFD)��Ŀ¼��Directory Entry(DE)��
˵����*******************************************************************************************************/
bool WriteTiff(const wchar_t * FileName,void *ptImage, size_t Height, size_t Width, size_t Byte)
{
	if (NULL==ptImage)
	{
		cout << "image open failed!" << endl;
		return false;
	}
	//CreateFile�������������CreateFileA��CreateFileW���ڴ���CreateFileWֱ��ע��
	HANDLE handle = CreateFileW(FileName, GENERIC_READ | GENERIC_WRITE, FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
	if (INVALID_HANDLE_VALUE==handle)
	{
		cout << "Create  FileName  Failed��" << endl; 
		return false;
	}
	//дIFH 8�ֽ�
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
	DWORD Length = 8+Byte*Height*Width; //��һ��IFD��λ��  --�����ļ���ָ��8+ptimage�Ĵ�С
	// Write the tiff file offset
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	if (4!=WriteSize)
	{
		cout << "Write the Offset Failed!" << endl;
		return false;
	}
	//д��Ӱ�������
	WriteFile(handle, ptImage, Byte*Width*Height, &WriteSize, NULL);
	//д��IFD�ļ�
	//DE����
	//���ߣ��ֱ��ʣ�����λ�ã��Ƿ�ѹ������ɫ���Լ�����Ϣ һ�����߸���Ϣ
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
	//д��Ӱ��Ŀ�
	tag = 0x0100;
	type = 0x0003;
	Length = 0x01;
	ValueOffset = Width;
	WriteFile(handle, &tag, 2, &WriteSize, NULL);
	WriteFile(handle, &type, 2, &WriteSize, NULL);
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	WriteFile(handle, &ValueOffset, 4, &WriteSize, NULL);
	//д��Ӱ��ĸ�
	tag = 0x0101;
	type = 0x0003;
	Length = 0x01;
	ValueOffset = Height;
	WriteFile(handle, &tag, 2, &WriteSize, NULL);
	WriteFile(handle, &type, 2, &WriteSize, NULL);
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	WriteFile(handle, &ValueOffset, 4, &WriteSize, NULL);
	//д��Ӱ��ķֱ���
	tag = 0x0102;
	type = 0x0003;
	Length = 0x01;       //�������ݸ����������ɫ��>2��Ϊ��ɫ��
	ValueOffset = Byte*8;//λ�ֵ=1Ϊ��ɫ��=4Ϊ16ɫ��=8Ϊ256ɫ��
	WriteFile(handle, &tag, 2, &WriteSize, NULL);
	WriteFile(handle, &type, 2, &WriteSize, NULL);
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	WriteFile(handle, &ValueOffset, 4, &WriteSize, NULL);
	//д���Ƿ�ѹ��(��һ�����岻��û��Ҳ���Բ�ѹ�����)
	tag = 0x0103;
	type = 0x0003;
	Length = 0x01;
	ValueOffset = 0x08;//��05��ʾѹ��������ѹ����
	WriteFile(handle, &tag, 2, &WriteSize, NULL);
	WriteFile(handle, &type, 2, &WriteSize, NULL);
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	WriteFile(handle, &ValueOffset, 4, &WriteSize, NULL);
	//д��Ӱ���Ƿ�ɫ��
	/*������һ����Ϣ�Ļ�������ļ�����ԭʼͼ��ķ�ɫ���*/
	tag = 0x0106;
	type = 0x0003;
	Length = 0x01;
	ValueOffset = 0x01;//��01��ʾ��ɫ�����򲻷�ɫ��
	WriteFile(handle, &tag, 2, &WriteSize, NULL);
	WriteFile(handle, &type, 2, &WriteSize, NULL);
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	WriteFile(handle, &ValueOffset, 4, &WriteSize, NULL);
	//д���������ڵ�λ��
	tag = 0x0111;
	type = 0x0004;
	Length = 0x01;
	ValueOffset = 0x08;
	WriteFile(handle, &tag, 2, &WriteSize, NULL);
	WriteFile(handle, &type, 2, &WriteSize, NULL);
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	WriteFile(handle, &ValueOffset, 4, &WriteSize, NULL);
	//д���Լ�����Ϣ
	tag = 0x013B;
	type = 0x0002;
	Length = 8;//��Ϣ�ֽ���+1
	//8  ǰ8���ֽ�;��*��*�ֽ��� ���صĴ�С;2  Directory Entry Countռ�õ�2�ֽ�;DE������7*12  ÿ��IFD��ռ���ֽ��� ;4  ��һ��tag��λ��Offset to next IFD;
	ValueOffset = 8+Height*Width*Byte+2+12*7+4;
	WriteFile(handle, &tag, 2, &WriteSize, NULL);
	WriteFile(handle, &type, 2, &WriteSize, NULL);
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	WriteFile(handle, &ValueOffset, 4, &WriteSize, NULL);
	//д���ļ�β
	Length = 0x00;
	WriteFile(handle, &Length, 4, &WriteSize, NULL);
	//д���Լ�����Ϣ
	char *text="YANGHUI";
	WriteFile(handle,text,8,&WriteSize,NULL);
	cout<<"������Ӱ�������ɣ�"<<endl;
	return true;
}

 



























