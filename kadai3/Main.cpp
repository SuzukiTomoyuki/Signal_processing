#include <vector>
#include <iostream>
#include <fstream>
#include <boost/array.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define PI 3.14159265358979323846
#define N 2

using namespace std;

//関数プロトタイプ宣言
void sin_init(float freq, float sps, float phi);
double sin_get();
void Smoothing();

boost::array<double, 1000> ar;
boost::array<double, N+1> ar_temp;
boost::array<double, 1000> ar_Smoothing;

double ar1,yz1,yz2;
double buf=0.0;
int count_smothing=0;

//初期化（freq=周波数　sps=サンプリング周波数 phi=位相）
void sin_init(float freq, float sps, float phi) {
	yz1 = sin(phi);
	yz2 = sin(phi+2*PI*freq/sps);
	ar1 = 2*cos(2*PI*freq/sps);

}

double sin_get() {
	double g;
	g   = yz1;
	yz1 = yz2;
	yz2 = ar1*yz1 - g;
	return g;
}

void Smoothing(){
	int i,j;
	double temp=0;

	for(i=0;i<1000;++i){
		for(j=-(N/2);j<N/2;++j){
			if((i+j)<0)temp+=ar[0];
			else
				if((i+j)>999)temp+=ar[999];
			else
				temp+=ar[i+j];
		}
		ar_Smoothing[i]=temp/(N+1);
		temp=0;
	}

}

int main()
{
	int i;
	sin_init(10, 1000, 0);

	clock_t start,end;
	start = clock();
	end = clock();

	std::ofstream ofs( "output.txt" );

	for(i=0;i<=999;++i){
		ar[i] = 5*sin_get()+(double) rand() / RAND_MAX;
	}

	Smoothing();

	for(i=0;i<=999;++i){
		ofs <<ar_Smoothing[i]<< std::endl;
		end = clock();
	}

	return 0;
}