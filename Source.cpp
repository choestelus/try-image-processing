#include <cstdio>
#include <climits>
#include <opencv\cv.h>
#include <opencv\highgui.h>
#include <opencv2\imgproc\imgproc.hpp>
#include <opencv2\video\video.hpp>
#include <opencv\highgui.h>
#include <iostream>
#include <vector>
#include <iterator>
#include <iomanip>
#include <algorithm>
#include <cstring>
#include "cvplot.h"
#include "statmath.h"
#define rowPtr(imagePtr, dataType, lineIndex) \
	(dataType *)(imagePtr->imageData + (lineIndex) * imagePtr->widthStep)

using namespace cv;
using namespace std;

void accessMat(Mat&, Mat&);

template <typename T>
int cumulative(T *in, int size);

template <class T>
int maxMat(T *in, int size);
template <class T>
T maxt(const T &a, const T &b)
{	return a>b?a:b;	}
int main(int argc, char** argv)
{	
	unsigned short max = 0;
	double maxmean = 0;
	Mat matimg;
	if(argc == 1)
		matimg = imread("crop_0094.tif",CV_LOAD_IMAGE_UNCHANGED);
	else
		matimg = imread(argv[1],CV_LOAD_IMAGE_UNCHANGED);
	if(!matimg.data)
	{
		cout<<"unable to load image\n";
		exit(0);
	}

	Mat gimg(matimg.rows, matimg.cols, CV_16U);
	//convert to grayscale
	accessMat(matimg, gimg);
	//comparison between grayscale and original
	namedWindow("image to graph");
	imshow("image to graph", matimg);
	namedWindow("grayscale");
	imshow("grayscale", gimg);

	// specify a line to plot
	int the_line = 0;
	int key = -1;
	while (the_line < gimg.rows)
	{
		unsigned short *pb = gimg.ptr<unsigned short>(the_line);
		vector<double> zpb;
		
		
		
		int width = gimg.cols;
		int heigth = gimg.rows;
		double mean = avrg(pb, width);
		double sd = stdv(pb, width, mean);
		int localmax = peakcount(pb, width,100);
		max = maxt<unsigned short>(maxMat<unsigned short>(pb, width),max);
		maxmean = maxt<double>(mean, maxmean);
		zpb.resize(width);


		for(unsigned int i=0; i<width; i++)
			zpb[i] = zscore(pb[i], mean, sd);

		cout<<"row  : ["<<the_line+1<<"/"<<heigth<<"]";
		cout<<" max its : "<<maxMat<unsigned short>(pb, width);
		cout<<" max OIS : "<<max;
		cout<<" mean : "<<fixed<<setprecision(2)<<mean;
		cout<<" maxmean : "<<fixed<<setprecision(2)<<maxmean;
		cout<<" s.d. : "<<fixed<<setprecision(2)<<sd;
		cout<<" zmin : "<<fixed<<setprecision(2)<<*min_element(zpb.begin(), zpb.end());
		cout<<" zmax : "<<fixed<<setprecision(2)<<*max_element(zpb.begin(), zpb.end());
		cout<<" pc : "<<localmax;
		cout<<endl;
		CvPlot::plot("GC", pb, width, 1, 128, 192, 128);
		//CvPlot::plot("memcpy", cpb, width, 1, 128, 192, 128);
		
		key = cvWaitKey(0);
		if (key == 32)
		{
			// plot the next line
			the_line++;
			// clear previous plots
			CvPlot::clear("GC");
			//CvPlot::clear("memcpy");
		}
		else
			break;
	}
	waitKey(0);
	return 0;
}

void accessMat(Mat &M, Mat &out)
{
	//unsigned short *p;
	unsigned short *it;
	unsigned short *sp[3];

	vector<Mat> splittedimg(3);
	split(M, splittedimg);
	for(int i=0; i<M.rows; i++)
	{
		it = out.ptr<unsigned short>(i);
		for(int pt=0; pt < 3; pt++)
			sp[pt] = splittedimg[pt].ptr<unsigned short>(i);
		for(int j=0; j<M.cols; j++)
		{
			it[j] = unsigned short(sp[0][j]*0.0722 + sp[1][j]*0.7152 + sp[2][j]*0.2126);
		}
	}
}

template<typename T>
int cumulative(T *in, int size)
{
	int sum = 0;
	for(int i=0;i < size; i++)
		sum+=in[i];
	return sum;
}

template <class T>
int maxMat(T *in, int size)
{
	int max = INT_MIN;
	for(int i=0; i<size; i++)
	{
		if(max < in[i])
			max = in[i];
	}
	return max;
}