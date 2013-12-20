#include <iostream>
#include <fstream>
#include <boost/array.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define PI 3.14159265358979323846	//円周率
#define N 1024		//要素数

using namespace std;

//関数プロトタイプ宣言
void sin_init(double freq, double sps, double phi);
double sin_get();
void dft();
void fft();

boost::array<double, N> ar;
boost::array<double, N> fft_re;
boost::array<double, N> fft_im;

double ar1,yz1,yz2;
double x[N], y[N], wx[N], wy[N]; 

//初期化（freq=周波数　sps=サンプリング周波数 phi=位相）
void sin_init(double freq, double sps, double phi) {
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

//FFTフーリエ変換
void fft(){ 

	std::ofstream ofs( "output2.txt" );

	int temp=N,m=0;
	double xx,xy; 
	int i,ii,il,j,jj,k,l,l2,nv2; 
	double z; 

	while(temp!=1){
		temp/=2;
		m++;
	} 

	nv2=N/2; 
	j=0;
	for(i=0;i<N-1;i++) { 
		if(i<j) { 
			z=x[j]; 
			x[j]=x[i]; 
			x[i]=z; 
			z=y[j]; 
			y[j]=y[i]; 
			y[i]=z; 
		} 
		k=nv2; 
		if(k<j+1) 
			do {
				j=j-k;
				k=k/2;
			} 
			while(k<j+1); 
			j=j+k;  
		} 
		for (i=1; i<=m; i++) { 
			l=1; 
			l=l<<i;
			l2=l/2;
			jj=0; 
			ii=N/l; 
			for (j=0; j<l2; j++) {
				for (k=j; k<N; k=k+l) { 
					il=k+l2;
					xx=x[il]*wx[jj]-y[il]*wy[jj];
					xy=y[il]*wx[jj]+x[il]*wy[jj];
					x[il]=x[k]-xx; 
					y[il]=y[k]-xy;
					x[k]=x[k]+xx;
					y[k]=y[k]+xy;
				}  
				jj=jj+ii; 
			} 
		} 
}

//DFTフーリエ変換
void dft(){
	std::ofstream ofs( "output.txt" );

	int i,j;
	double im=0,re=0;

	for(i=0;i<N;++i){
		re=im=0;
		for(j=0;j<N;++j){
			re+=ar[j]*cos(2*PI*i*j/N);
			im+=-(ar[j]*sin(2*PI*i*j/N));
		}
		//実部、虚部の順にファイルに書き出す
		ofs <<re<<'\t'<<im<< std::endl;
	}
}

int main(){

	std::ofstream ofs( "output2.txt" );
	int i;

	sin_init(10, 1000, 0);

	clock_t start,end;
	start = clock();
	end = clock();

	for(i=0;i<N;++i){
		//ar[i]=5*sin_get();
		ar[i]=5*sin_get();//+(double) rand() / RAND_MAX;
	}

	dft();

	for (i=0; i<N; i++){ 
		wx[i]=cos(2*PI*i/N); 
		wy[i]=-sin(2*PI*i/N);
		x[i]=5*sin_get();
		y[i]=0; 
	} 
	x[1]=1;




	fft();
	for (i=0;i<N; i++){ 
		//printf("%f\t%f\n",x[i],y[i]);
		//実部、虚部の順にファイルに書き出す
		ofs <<x[i]<<'\t'<<-y[i]<< std::endl;
		} 

	return 0;
}