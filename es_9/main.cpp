#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "GA.h"
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
	int nc=32;
	int np=900;
	
	//CIRCONFERENZA
	GA ga1(nc, np);
	ga1.GenerateCities(0);
	ga1.GeneratePopulation();
	ofstream x1;
	x1.open("L_best_1");
	ofstream y1;
	y1.open("L_ave_1");
	ofstream z1;
	z1.open("L_err_1");
	ofstream b1;
	b1.open("best_1");
	ofstream c1;
	c1.open("city_1");
	for (int i=0; i<1000; i++){
		ga1.NewGeneration();
		ga1.NLBest(int(np/2));
		double L=media(ga1.nL_best, int(np/2));
		double L_err=sqrt(varianza(ga1.nL_best, int(np/2)));
		double L_best=ga1.nL_best[0];	
		if (x1.is_open()){x1<<L_best<<endl;}
		else cerr<<"problema: non si è aperto il file "<<endl;
		if (y1.is_open()){y1<<L<<endl;}
		else cerr<<"problema: non si è aperto il file "<<endl;
		if (z1.is_open()){z1<<L_err<<endl;}
		else cerr<<"problema: non si è aperto il file "<<endl;
		//cout<<"circonferenza: "<<L_best<<"  "<<L<<" +- "<<L_err<<endl;
	}
	for (int i=0; i<nc; i++){
		if (b1.is_open()){b1<<ga1.population[ga1.iBest()][i]<<endl;}
		else cerr<<"problema: non si è aperto il file "<<endl;
		if (c1.is_open()){c1<<ga1.cities[i][0]<<" "<<ga1.cities[i][1]<<endl;}
		else cerr<<"problema: non si è aperto il file "<<endl;
	}
	b1<<0<<endl;
	x1.close();
	y1.close();
	z1.close();
	b1.close();
	c1.close();
	
	//QUADRATO
	GA ga2(nc, np);
	ga2.GenerateCities(1);
	ga2.GeneratePopulation();
	ofstream x2;
	x2.open("L_best_2");
	ofstream y2;
	y2.open("L_ave_2");
	ofstream z2;
	z2.open("L_err_2");
	ofstream b2;
	b2.open("best_2");
	ofstream c2;
	c2.open("city_2");
	for (int i=0; i<1000; i++){
		ga2.NewGeneration();
		ga2.NLBest(int(np/2));
		double L=media(ga2.nL_best, int(np/2));
		double L_err=sqrt(varianza(ga2.nL_best, int(np/2)));
		double L_best=ga2.nL_best[0];	
		if (x2.is_open()){x2<<L_best<<endl;}
		else cerr<<"problema: non si è aperto il file "<<endl;
		if (y2.is_open()){y2<<L<<endl;}
		else cerr<<"problema: non si è aperto il file "<<endl;
		if (z2.is_open()){z2<<L_err<<endl;}
		else cerr<<"problema: non si è aperto il file "<<endl;
		cout<<"rettangolo: "<<L_best<<"  "<<L<<" +- "<<L_err<<endl;
	}
	for (int i=0; i<nc; i++){
		if (b2.is_open()){b2<<ga2.population[ga2.iBest()][i]<<endl;}
		else cerr<<"problema: non si è aperto il file "<<endl;
		if (c2.is_open()){c2<<ga2.cities[i][0]<<" "<<ga2.cities[i][1]<<endl;}
		else cerr<<"problema: non si è aperto il file "<<endl;
	}
	b2<<0<<endl;
	x2.close();
	y2.close();
	z2.close();
	b2.close();
	c2.close();
	
	
	return 0;
}
