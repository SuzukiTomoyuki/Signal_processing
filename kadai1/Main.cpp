#include <iostream>
#include <fstream>
#include <boost/array.hpp>
#include <stdio.h>
#include <math.h>
#include <time.h>
#define PI 3.14159265358979323846

using namespace std;

float ar1,yz1,yz2;

//初期化（freq=周波数　sps=サンプリング周波数 phi=位相）
void sin_init(float freq, float sps, float phi) {
	yz1 = sin(phi);
	yz2 = sin(phi+2*PI*freq/sps);
	ar1 = 2*cos(2*PI*freq/sps);

}

float sin_get() {
	float g;
	g   = yz1;
	yz1 = yz2;
	yz2 = ar1*yz1 - g;
	return g;
}

int main()
{
	int i;
	double x=127/5;
	boost::array<double, 1000> ar;

	sin_init(10, 1000, 0);

	clock_t start,end;
	start = clock();
	end = clock();

	std::ofstream ofs( "output.txt" );

	for(i=0;i<=999;++i){
		ar[i]=5*sin_get();
		ofs <<ar[i]<<'\t'<<(int)(ar[i]*x)<< std::endl;
		end = clock();
	}

	return 0;
}