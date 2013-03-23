#pragma once
#ifndef _statmath_h
#define _statmath_h
#include <vector>
#include <cmath>
#include <algorithm>
const double PI = 4*atan(1);

template <class T>
T maxt(const T &a, const T &b)
{	return a>b?a:b;	}
double avrg(std::vector<unsigned short> v)
{      double sum=0;
       for(unsigned int i=0;i<v.size();i++)
               sum+=v[i];
       return sum/v.size();
}
double avrg(std::vector<double> v)
{
	double sum = 0;
    for(unsigned int i=0;i<v.size();i++)
		sum+=v[i];
	return sum/v.size();
}
double avrg(unsigned short *p, int num)
{
	double sum = 0;
	for(int i=0; i<num; i++)
		sum += p[i];
	return double(sum)/double(num);
}
double avrg(unsigned short *p, int init, int end)
{
	double sum = 0;
	for(int i=init; i<end; i++)
		sum += p[i];
	return double(sum)/double(end - init);
}
double gavrg(unsigned short *p, int num)
{
	return -1;
}
double stdv(std::vector<unsigned short> v, double ave)
{
    double E=0;
    // Quick Question - Can vector::size() return 0?
    double inverse = 1.0 / static_cast<double>(v.size());
    for(unsigned int i=0;i<v.size();i++)
    {
        E += pow(static_cast<double>(v[i]) - ave, 2);
    }
    return sqrt(inverse * E);
}
double stdv(unsigned short *p, int num, double mean)
{
	double E = 0;
	double inverse = 1.0 / static_cast<double>(num);
	for(int i = 0; i<num; i++)
		E += pow(static_cast<double>(p[i]) - mean, 2);
	return sqrt(inverse * E);
}
double zscore(unsigned short i, double mean, double stdv)
{
	return (double(i)-mean)/stdv;
}
double gfilter(double input, double sd)
{
	cout<<-(input*input)/(2*sd*sd)<<endl;
	cout<<exp(-(input*input)/(2*sd*sd))<<endl;
	cout<<(sqrt(2*PI)*sd)<<endl;
	return exp(-(input*input)/(2*sd*sd))/(sqrt(2*PI)*sd);
}
void gaussian1d(unsigned short *p, int num, double sd)
{
	for(int i=0; i<num; i++)
		p[i] = unsigned short(gfilter(double(p[i]),sd));
}
void gaussian1d(double *p, int num, double sd)
{
	for(int i=0; i<num; i++)
	{
		if(i == 100)
			cout<<"gaussian :"<<p[i]<<" "<<gfilter(double(p[i]),sd)<<endl;
		p[i] = gfilter(double(p[i]),sd);
	}
}
int peakcount(unsigned short *in, int num)
{
	int count = 0;
	for(int i=1;i<num-1;i++)
	{
		if(in[i-1] < in[i] && in[i] > in[i+1])
			count++;
	}
	return count;
}
int peakcount(unsigned short *in,const int num, int framesize, unsigned short *out)
{
	int count = 0;
	int wframsize = framesize;
	int worked = 0;
	while(worked < num)
	{
		double mean = avrg(in, worked, worked + wframsize);
		//cout<<"fmean "<<mean<<" ";
		for(int j=worked; j<worked+wframsize; j++)
		{
			if( in[j] > in[j-1] && in[j] > in[j+1] && in[j] >= mean)
			{
				out[j] = 40000;
				count++;
			}
			else
				out[j] = 20000;
		}
		if(worked + wframsize > num)
			wframsize = num - worked;
		worked += wframsize;
		//cout<<"wked "<<worked<<"num "<<num<<" \n";
	}
	return count;
}
int peakcount(unsigned short *in,const int num, int framesize)
{
	int count = 0;
	int wframsize = framesize;
	int worked = 0;
	while(worked < num)
	{
		double mean = avrg(in, worked, worked + wframsize);
		//cout<<"fmean "<<mean<<" ";
		for(int j=worked; j<worked+wframsize; j++)
		{
			if( in[j] > in[j-1] && in[j] > in[j+1] && in[j] >= mean)
				count++;
		}
		if(worked + wframsize > num)
			wframsize = num - worked;
		worked += wframsize;
		//cout<<"wked "<<worked<<"num "<<num<<" \n";
	}
	return count;
}
int peakcount(unsigned short *in,const int num, double sd)
{
	int count = 0;
	double mean = avrg(in, num);
	for(int i=0; i<num; i++)
	{
		if(in[i] > mean+sd)
			count++;
	}
	return count;
}
int peakcount(std::vector<double> zvec, bool *out)
{
	int count = 0;
	double zsum = 0, zmean = 0;
	int zc = 0;
	for(unsigned int i=0;i<zvec.size(); i++)
	{
		if(zvec[i] > 0)
		{
			zsum += zvec[i];
			zc++;
		}
	}
	zmean = zsum/zc;
	//cout<<"zm"<<zmean<<"zc"<<zc<<" ";
	for(unsigned int i=0;i<zvec.size(); i++)
	{
		if(zvec[i] > zmean)
		{
			count++;
			out[i] = true;
		}
	}
	return count;
}
void designal(double *out, unsigned short *in, double *insignal, int num)
{
	for(int i=0; i<num; i++)
		out[i] = double(in[i]) - insignal[i];
}
void simpmovavg(unsigned short *in, double *out, int num, int n)
{
	//initializing out of range data
	for(int i=0; i<n; i++)
		out[i] = 0;
	for(int i=0; i<n; i++)
	{
		if(i == 0)
			out[i] = in[i];
		else
		{
			int tempsum = 0;
			for(int j = i; j<n;j++)
				tempsum += in[j];
			out[i] = double(tempsum)/double(i+1);
		}
	}
	//end of out of range case

	double avgsum = 0;
	int curr = 0;

	for(int i=curr; i<curr+n; i++)
		avgsum += in[i];
	avgsum = avgsum/n;
	out[n-1] = avgsum;
	curr = n;

	while(curr < num)
	{
		out[curr] = out[curr-1] - double(in[curr-1])/double(n) + double(in[curr])/double(n);
		curr++;
	}
}
void simpmovmed(unsigned short *in, double *out, int num, int n)
{
	int curr = 0;
	if(n%2 == 0)
		n++;
	std::vector<unsigned short> temp;
	temp.resize(n);

	while(curr < num)
	{
		//left out of range case
		if( curr < n/2 )
		{
			//left out of range element
			for(int i=0; i < n/2 - curr; i++)
				temp[i] = in[0];
			//mid element and beyond; n/2 - curr is difference between two iterator
			for(int i = n/2 - curr; i < n - n/2 + curr; i++)
				temp[i] = in[i + n/2 - curr];

			std::sort(temp.begin(),temp.end());
			out[curr] = temp[n/2];
			curr++;
		}
		//right out of range case
		else if(curr + n/2  >= num)
		{
			for(int i = n - (curr+(n/2)-num-1);i < n;i++)
				temp[i] = in[num-1];
			for(int i = 0; i < n - (curr+(n/2)-num-1); i++)
				temp[i] = in[i + curr - n/2];
			curr++;
		}
		//general case
		else
		{
			for(int i=curr-n/2; i<1+curr+n/2; i++)
				temp[i] = in[i];
			std::sort(temp.begin(), temp.end());
			out[curr] = temp[n/2];
			curr++;
		}
	}

}
#endif