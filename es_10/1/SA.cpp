#include "SA.h"
#define _USE_MATH_DEFINES 

using namespace std;

SA::SA(int Nc, double Tstart, double deltaT, int Tstep, int step){			//inizializzazione
	_Nc = Nc;
	_Tstart = Tstart;
	_deltaT=deltaT;
	_Tstep = Tstep;
	_step = step;
	
	cities.resize(_Nc);
	for(int i = 0; i < _Nc; i++) {cities[i].resize(2);}
	cammino.resize(_Nc);
	newcammino.resize(_Nc);
	bestcammino.resize(_Nc);
	
	
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){Primes >> p1 >> p2 ;} 
	else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();
	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} 
	else cerr << "PROBLEM: Unable to open seed.in" << endl;
}

SA::~SA(){}				//rimozione

void SA::GenerateCities(int mode){			//generazione città nelle due modalità
	if(mode==0){
		for(int i = 0; i < _Nc; i++){
			double teta = rnd.Rannyu(0., 2.*M_PI);
			double xc = cos(teta);
			double yc = sin(teta);
			cities[i][0]=xc;
			cities[i][1]=yc;
		}
	}
	else if(mode==1){
		for(int i = 0; i < _Nc; i++){
			cities[i][0]=rnd.Rannyu(-1., 1.);
			cities[i][1]=rnd.Rannyu(-1., 1.);
		}
	}
}

double SA :: Distanza(vector<int> cammino){
	double dist=0;						//calcolo distanze salvate in dist
	for (int i=0; i<_Nc; i++){
		int c1=cammino[i];
		int c2=cammino[(i+1)%_Nc];
		double d=sqrt(pow(cities[c1][0]-cities[c2][0],2)+pow(cities[c1][1]-cities[c2][1],2));
		dist=dist+d;
	}
	return dist;
}

void SA::GenerateCammino(){			//generazione popolazione
	cammino[0]=0;
	bestcammino[0]=0;
	for(int i=1; i<_Nc; i++){
		int r;
		int control;
		do{
			control=0;
			r= int(_Nc*rnd.Rannyu());
			for(int j=0; j<i; j++){
				if(cammino[j]==r){control=1;}
			}
		}
		while(control==1);
		cammino[i]=r;
		bestcammino[i]=r;
		
	}
	rnd.SaveSeed();
	_dist=Distanza(cammino);
	_bestdist=Distanza(bestcammino);
}

void SA :: Mutation(int type=0){
	if (type==0){					//permutazione di due città
		int a=int((_Nc-1)*rnd.Rannyu())+1;
		int b=int((_Nc-1)*rnd.Rannyu())+1;
		int c=newcammino[a];
		newcammino[a]=newcammino[b];
		newcammino[b]=c;
	}
	else if (type==1){			//shift di ncit città di npos posizioni a partire dalla città ic
		vector<int> temp(_Nc-1);
		for (int j=0; j<(_Nc-1); j++){temp[j]=newcammino[j+1];}
		int ncit=(_Nc-1)*rnd.Rannyu();	
		int npos=(_Nc-1-ncit)*rnd.Rannyu();
		int ic=(_Nc-1)*rnd.Rannyu();
		vector<int> t_cit(ncit);
		vector<int> t_pos(npos);
		for (int j=0; j<ncit; j++){t_cit[j]=temp[(ic+j)%(_Nc-1)];}
		for (int j=0; j<npos; j++){t_pos[j]=temp[(ic+ncit+j)%(_Nc-1)];}
		for (int j=0; j<npos; j++){temp[(ic+j)%(_Nc-1)]=t_pos[j];}
		for (int j=0; j<ncit; j++){temp[(ic+npos+j)%(_Nc-1)]=t_cit[j];}
		for (int j=0; j<(_Nc-1); j++){newcammino[j+1]=temp[j];}
	}
	if (type==2){						//scambio di ncit adiacenti
		vector<int> temp(_Nc-1);
		for (int j=0; j<(_Nc-1); j++){temp[j]=newcammino[j+1];}
		int ncit=((_Nc-1)/2)*rnd.Rannyu();	
		int ic=(_Nc-1)*rnd.Rannyu();
		vector<int> t_cit(ncit);
		for (int j=0; j<ncit; j++){t_cit[j]=temp[(ic+j)%(_Nc-1)];}
		for (int j=0; j<ncit; j++){temp[(ic+j)%(_Nc-1)]=temp[(ic+ncit+j)%(_Nc-1)];}
		for (int j=0; j<ncit; j++){temp[(ic+ncit+j)%(_Nc-1)]=t_cit[j];}
		for (int j=0; j<(_Nc-1); j++){newcammino[j+1]=temp[j];}
	}
	if (type==3){						//inversione ordine di ncit adiacenti
		vector<int> temp(_Nc-1);
		for (int j=0; j<(_Nc-1); j++){temp[j]=newcammino[j+1];}
		int ncit=(_Nc-1)*rnd.Rannyu();	
		int ic=(_Nc-1)*rnd.Rannyu();
		vector<int> t_cit(ncit);
		for (int j=0; j<ncit; j++){t_cit[j]=temp[(ic+j)%(_Nc-1)];}
		for (int j=0; j<ncit; j++){temp[(ic+j)%(_Nc-1)]=t_cit[(ncit-1-j)];}
		for (int j=0; j<(_Nc-1); j++){newcammino[j+1]=temp[j];}
	}

}

void SA :: EvolutionT(double T){
	for (int i=0; i<_Nc; i++){newcammino[i]=cammino[i];}
	Mutation(0);
	Mutation(1);
	Mutation(2);
	Mutation(3);
	_dist=Distanza(cammino);
	_newdist=Distanza(newcammino);
	double p=exp(-T*(_newdist-_dist));
	if(_newdist<_dist){p=1;}
	if(rnd.Rannyu()<p){
		for(int i=0; i<_Nc; i++){cammino[i]=newcammino[i];}
		_dist=_newdist;
	}
	if(_newdist<_bestdist){
		for(int i=0; i<_Nc; i++){bestcammino[i]=newcammino[i];}
		_bestdist=_newdist;
	}
}

void SA :: Evolution(const char *nome){
	ofstream x;
	x.open(nome);
	double T=_Tstart;
	for (int j=0; j<_Tstep; j++){
		cout<<"current temperature (che in realtà è beta): "<<T<<endl;
		for (int i=0; i<_step; i++){
			EvolutionT(T);
		}
		x<<_bestdist<<" "<<T<<endl;
		T=T+_deltaT;
	}
}
















