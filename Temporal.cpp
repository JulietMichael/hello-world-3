//This code outputs data for the last three periods for a single specified frequency and field amplitude.
//Initial polarization and the first derivitive of the initial polarization is set to zero.

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
using namespace std;

//This aides in assigning output file names based on simulation parameters
template <typename T>
string NumberToString(T pNumber)
{
 ostringstream oOStrStream;
 oOStrStream << pNumber;
 return oOStrStream.str();
}

//pi is meant to represent the circumferene of a circule to its diameter and the rest of these values correspond
//to the coefficients of the 2nd order ODE solved by the code.
const double pi=3.14159265359;
double C1=3000;
double C2=15.24830262;
double C3=-4.71554993;
double C4=71.15337659;

//these are the functions corresponding to the two 1st order ODEs that are used to 
//break up the 2nd order ODE that this equation solves.
double u1(double u2);
double u2(double u1, double u2, double t, double E0, double f);

//number of data points per period
const int x=10000000;
	
//number of data points to skip per period when outpurring data
const int skips=1000;

//number of data points retained for outputputting per frequency per field amplitude
const int dataPts=x/skips;

//number of different frequencies to simulate per amplitude
const int trajectoriesF=1;

//number of different field amplitudes to simulate per frequency
const int trajectoriesE=1;

//thse contribute to file output file names
string extension=".txt";
string preextensionE="E";
string preextensionP="P";

int main()
{
	//variable that will be used by the Runga-Kutta algorithm
	double k11;
	double k12; 
	double k21; 
	double k22; 
	double k31; 
	double k32; 
	double k41; 
	double k42; 
	
	ofstream Es ;   
	ofstream Ps;
	ofstream info("f.txt");
	ofstream times("time.txt");
	
	//amplitude of E field
	double E0=100;

	for(int q=0;q<trajectoriesE;q++)
	{
		//multiplies the fundamental frequency
		double multiplier=150;

		string filename;
		filename=NumberToString(E0);		    
		Es.open((preextensionE+filename+extension).c_str());	
		Ps.open((preextensionP+filename+extension).c_str());

		for(int p=0;p<trajectoriesF;p++)
		{			
			//fundamental frequency (angular)
			double f=multiplier*(2*pi/6400);

			double period=6400/multiplier;

			double m=period/x;
			
			//initial conditions
			double w1=0;
			double w2=0;
			double t=0;
				
			//number of periods to examine
			int periods=6+multiplier;
			
			//max number of iterations
			int max=periods*x;
			
			info << multiplier << endl;
			
			int index=0;
			
			for(int i=0;i<=max;i++)
			{
				if(i==0)
				{
					w1=0;
					w2=0;
				}
				
				else
				{
				
					k11=m*u1(w2);
					k12=m*u2(w1, w2, t, E0, f);

					k21=m*u1(w2+.5*k12);
					k22=m*u2(w1+.5*k11, w2+.5*k12, t+m/2, E0, f);
					
					k31=m*u1(w2+.5*k22);
					k32=m*u2(w1+.5*k21, w2+.5*k22, t+m/2, E0, f);
				
					k41=m*u1(w2+k32);
					k42=m*u2(w1+k31, w2+k32, t+m, E0, f);
				
					w1=w1+(k11+2*k21+2*k31+k41)/6;
					w2=w2+(k12+2*k22+2*k32+k42)/6;
				
				}
				t=i*m;
				
				if((m*i)>(m*(max+1)-3*6400/multiplier) && i%skips==0)
				{	
					times << t << endl;
					Es << E0*sin(f*t) << endl;
					Ps << w1 << endl;
					index++;
				}	
			}

			multiplier++;
		}
		
		if(q==0)
		{
			info.close();
			times.close();
		}

		Es.close();
		Ps.close();		
		
		E0=E0+10;
	}
	
	return 0;
}

double u1(double u2)
{
	return u2;
}

double u2(double u1, double u2, double t, double E0, double f)
{
	double response;
	response=-C1*u2+C2*u1+E0*sin(f*t)-C3*pow(u1,3)-C4*pow(u1,5);
	return response;
}