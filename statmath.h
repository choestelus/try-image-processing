#pragma once
#ifndef _statmath_h
#define _statmath_h
#include <vector>
#include <cmath>
const double PI = 4*atan(1);
double avrg(std::vector<unsigned short> v)
{      int sum=0;
       for(unsigned int i=0;i<v.size();i++)
               sum+=v[i];
       return sum/v.size();
}
double avrg(unsigned short *p, int num)
{
	int sum = 0;
	for(unsigned int i=0; i<num; i++)
		sum += p[i];
	return double(sum)/double(num);
}
double avrg(unsigned short *p, int init, int end)
{
	int sum = 0;
	for(unsigned int i=init; i<end; i++)
		sum += p[i];
	return double(sum)/double(end - init);
}
double gavrg(unsigned short *p, int num)
{

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
	for(unsigned int i = 0; i<num; i++)
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
int peakcount(unsigned short *in,const int num, int framesize)
{
	//bool *printed = new bool[num];
	bool *cted = new bool[num];
	for(int i=0; i<num;i++)
	{
		//printed[i] = false;
		cted[i] = false;
	}
	int count = 0;
	//cout<<"[ ";
	for(int i=0; i<num-framesize; i+=1)
	{
		double mean = avrg(in, i, i+framesize);
		for(int j = i; j<i+framesize; j++)
		{
			if( in[j] > mean && in[j] > in[j-1] && in[j] > in[j+1] && !cted[j])
			{
				count++;
				cted[j] = true;
			}
		}
	}
	//cout<<" ]"<<endl;
	//delete[] printed;
	delete[] cted;
	return count;
}

/*double spectralflat(unsigned short *in, int num)
{
	double mean = avrg(in,num);
}*/
#endif