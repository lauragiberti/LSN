#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "SA.h"
#include <vector>
#define _USE_MATH_DEFINES 
using namespace std;

//media delle prime N componenti di un vettore v
double media (vector<double>v, int N) {
	double sum=0;
	for (int i=0; i<N; i++){
		sum=sum+v[i];
	}
	return sum/N;
}

//varianza delle prime N componenti di un vettore v
double varianza (vector<double>v, int N){
	double v1=0;
	double v2=0;
	for (int i=0; i<N; i++){
		v1+=v[i];
		v2+=v[i]*v[i];
	}
	v1=v1/N;
	v2=v2/N;
	return (v2-v1*v1);
}

int main(){
	int Nc=32;
	double Tstart=0;
	double deltaT=0.5;
	int Tstep=200;
	int step=10000;
	
	//CIRCONFERENZA
	SA sa (Nc, Tstart, deltaT, Tstep, step);	
	sa.GenerateCities(0);	
	sa.GenerateCammino();	
	sa.Evolution("L_best_1");
	
	ofstream x;
	
	x.open("city_1");
	for (int i=0; i<Nc; i++){x<<sa.cities[i][0]<<" "<<sa.cities[i][1]<<endl;}
	x.close();
	
	x.open("bestcammino_1");
	for (int i=0; i<Nc; i++){x<<sa.bestcammino[i]<<endl;}
	x<<0<<endl;
	x.close();
	
	
	//RETTANGOLO
	SA saa (Nc, Tstart, deltaT, Tstep, step);	
	saa.GenerateCities(1);	
	saa.GenerateCammino();	
	saa.Evolution("L_best_2");
	
	x.open("city_2");
	for (int i=0; i<Nc; i++){x<<saa.cities[i][0]<<" "<<saa.cities[i][1]<<endl;}
	x.close();
	
	x.open("bestcammino_2");
	for (int i=0; i<Nc; i++){x<<saa.bestcammino[i]<<endl;}
	x<<0<<endl;
	x.close();
	
	
	return 0;
}
