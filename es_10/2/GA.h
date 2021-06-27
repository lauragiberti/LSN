#ifndef __GA__
#define __GA__
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include "random.h"

//Genetic Algorithm class
using namespace std;

class GA{
	public:
	//protected:
	int Nc;
	int Np;
	Random rnd;
	vector<vector<double> > cities;
	vector<vector<int> > population;
	vector<vector<int>> newpopulation;
	vector<double> dist;
	vector<int> idist;
	vector<double> nL_best;
	
	//public:
	GA(int, int, int);
	~GA();
	//int GetNc(){return Nc;};
	//int GetNp(){return Np;};
	void GenerateCities(int);
	void GeneratePopulation();
	void Distanze();
	void Order();
	int Selection();
	void Crossover();
	void Mutation(int);
	void NewGeneration();
	void NLBest(int);
	int iBest();
};

#endif
