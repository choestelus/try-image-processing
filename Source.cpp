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
#define mframesize 25
using namespace cv;
using namespace std;
using std::cout;

void accessMat(Mat&, Mat&);

template <typename T>
int cumulative(T *in, int size);

// maximum value in 1D Matrix
template <class T>
int maxMat(T *in, int size);


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
	//namedWindow("image to graph");
	//imshow("image to graph", matimg);
	//namedWindow("grayscale");
	//imshow("grayscale", gimg);

	// specify a line to plot
	int the_line = 0;
	int key = -1;
	
	int maxrow = 0;
	double sd = 0;
	bool *peakmap = nullptr;
	bool pmallocated = false;
	bool movavgalloc = false;
	double *movavg = nullptr;
	double *noise = nullptr;
	while (the_line < gimg.rows)
	{
		unsigned short *pb = gimg.ptr<unsigned short>(the_line);
		vector<double> zpb;

		int width = gimg.cols;
		int heigth = gimg.rows;
		double mean = avrg(pb, width);
		double sd = stdv(pb, width, mean);

		if(!pmallocated)
		{
			peakmap = new bool[width];
			pmallocated = true;
		}
		if(!movavgalloc)
		{
			movavg = new double[width];
			noise = new double[width];
			movavgalloc = true;
		}
		//simpmovavg(pb, movavg, width, mframesize);

		for(int i=0; i<width; i++)
			peakmap[i] = false;
		max = maxt<unsigned short>(maxMat<unsigned short>(pb, width),max);
		if(mean > maxmean)
			maxrow = the_line;
		maxmean = maxt<double>(mean, maxmean);
		zpb.resize(width);

		for(int i=0; i<width; i++)
			zpb[i] = zscore(pb[i], mean, sd);
		
		/*cout<<"row  : ["<<the_line+1<<"/"<<heigth<<"]";
		cout<<" max its : "<<maxMat<unsigned short>(pb, width);
		cout<<" max OIS : "<<max;
		cout<<" mean : "<<fixed<<setprecision(2)<<mean;
		cout<<" maxmean : "<<fixed<<setprecision(2)<<maxmean;
		cout<<" s.d. : "<<fixed<<setprecision(2)<<sd;
		cout<<" zmin : "<<fixed<<setprecision(2)<<*min_element(zpb.begin(), zpb.end());
		cout<<" zmax : "<<fixed<<setprecision(2)<<*max_element(zpb.begin(), zpb.end());
		
		cout<<endl;
		CvPlot::plot("GC", pb, width, 1, 128, 192, 128);
		CvPlot::label("input");
		CvPlot::plot("GC", movavg, width, 1, 255, 0, 0);
		CvPlot::label("smovavg");
		key = cvWaitKey(0);
		if (key == 32)
		{
			// plot the next line
			the_line++;
			// clear previous plots
			CvPlot::clear("GC");
			//CvPlot::clear("peak");
		}
		else
			break;*/
		//in case of 1 line output
		the_line++;
	}

	{
		unsigned short *pb = gimg.ptr<unsigned short>(maxrow);
		vector<double> zpb;

		int width = gimg.cols;
		int heigth = gimg.rows;
		double mean = avrg(pb, width);
		double sd = stdv(pb, width, mean);

		zpb.resize(width);

		for(int i=0; i<width; i++)
			zpb[i] = zscore(pb[i], mean, sd);

		simpmovmed(pb, movavg, width, mframesize);
		designal(noise, pb, movavg, width);
		cout<<"test";
		cout<<"r"<<maxrow+1;
		cout<<" : mits "<<maxMat<unsigned short>(pb, width);
		cout<<" mmean "<<fixed<<setprecision(2)<<maxmean;
		cout<<" sd "<<fixed<<setprecision(2)<<sd;
		cout<<" z "<<fixed<<setprecision(2)<<*min_element(zpb.begin(), zpb.end());
		cout<<" ,";
		cout<<" "<<fixed<<setprecision(2)<<*max_element(zpb.begin(), zpb.end());
		cout<<endl;
		CvPlot::plot("GC", pb, width, 1, 128, 192, 128);
		CvPlot::label("input");
		CvPlot::plot("GC", movavg, width, 1, 255, 0, 0);
		CvPlot::label("signal");
		CvPlot::plot("Noise", noise, width, 1, 0, 0, 255);
	}
	if(movavgalloc)
	{
		delete[] movavg;
		delete[] noise;
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
// maximum value in 1D Matrix
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