#ifndef __SA__
#define __SA__
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

class SA{
	public:
	//protected:
	int _Nc;
	double _Tstart, _deltaT;
	int _Tstep, _step;
	Random rnd;
	double _dist, _newdist, _bestdist;
	vector<vector<double> > cities;
	vector<int> cammino;
	vector<int> newcammino;
	vector<int> bestcammino;
		
	//public:
	SA(int, double, double, int, int);
	~SA();
	//int GetNc(){return Nc;};
	//int GetNp(){return Np;};
	void GenerateCities(int);
	double Distanza(vector<int>);
	void GenerateCammino();
	void Mutation(int);
	void EvolutionT(double);
	void Evolution(const char *);
};

#endif
