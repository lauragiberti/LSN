#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "GA.h"
#include <vector>
#include "mpi.h"
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

int main(int argc, char* argv[]){

	int size, rank;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	cout<<"nodo: "<<rank<<" di: "<<size<<endl;
	
	int nc=32;
	int np=900;
	int Nmigr=10;

	GA ga(nc, np, rank);
	ga.GenerateCities(1);
	ga.GeneratePopulation();
	
	
	for (int i=0; i<50; i++){
		double L, L_best;
		vector <double> best_path0;
		best_path0.resize(nc);
		vector <double> best_path1;
		best_path1.resize(nc);
		vector <double> best_path2;
		best_path2.resize(nc);
		vector <double> best_path3;
		best_path3.resize(nc);
		
		for (int j=0; j<Nmigr; j++){
			ga.NewGeneration();
			ga.NLBest(10);
			L=media(ga.nL_best, 10);
			L_best=ga.nL_best[0];	
		}
		cout<<"generazione: "<<i<<" rank: "<<rank<<" -> "<<L_best<<"  "<<L<<endl;
		

		if (rank==0){
			for (int k=0; k<nc; k++){best_path0[k]=ga.population[ga.iBest()][k];}
		}
		if (rank==1){
			for (int k=0; k<nc; k++){best_path1[k]=ga.population[ga.iBest()][k];}
		}
		if (rank==2){
			for (int k=0; k<nc; k++){best_path2[k]=ga.population[ga.iBest()][k];}
		}
		if (rank==3){
			for (int k=0; k<nc; k++){best_path3[k]=ga.population[ga.iBest()][k];}
		}
		
		MPI_Bcast(&best_path0, nc, MPI_INTEGER, 0, MPI_COMM_WORLD);
		MPI_Bcast(&best_path1, nc, MPI_INTEGER, 1, MPI_COMM_WORLD);
		MPI_Bcast(&best_path2, nc, MPI_INTEGER, 2, MPI_COMM_WORLD);
		MPI_Bcast(&best_path3, nc, MPI_INTEGER, 3, MPI_COMM_WORLD);
		
	}
	
	
	return 0;
}
