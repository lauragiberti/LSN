#include "GA.h"
#define _USE_MATH_DEFINES 

using namespace std;

GA::GA(int nc, int np, int rank){			//inizializzazione
	Nc = nc;
	Np = np;
	
	population.resize(Np);
	for(int i = 0; i < Np; i++) {population[i].resize(Nc);}
	
	newpopulation.resize(Np);
	for(int i = 0; i < Np; i++) {newpopulation[i].resize(Nc);}
	
	cities.resize(Nc);
	for(int i = 0; i < Nc; i++) {cities[i].resize(2);}
	
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
	
	//if (rank==1){rnd.Rannyu();}
	//if (rank==2){rnd.Rannyu(); rnd.Rannyu();}
	//if (rank==3){rnd.Rannyu(); rnd.Rannyu(); rnd.Rannyu();}
}

GA::~GA(){}				//rimozione

void GA::GenerateCities(int mode){			//generazione città nelle due modalità
	if(mode==0){
		for(int i = 0; i < Nc; i++){
			double teta = rnd.Rannyu(0., 2.*M_PI);
			double xc = cos(teta);
			double yc = sin(teta);
			cities[i][0]=xc;
			cities[i][1]=yc;
		}
	}
	else if(mode==1){
		for(int i = 0; i < Nc; i++){
			cities[i][0]=rnd.Rannyu(-1., 1.);
			cities[i][1]=rnd.Rannyu(-1., 1.);
		}
	}
}

void GA::GeneratePopulation(){			//generazione popolazione
	for(int i=0; i <Np; i++){
		population[i][0]=0;
		for(int j=1; j<Nc; j++){
			int r;
			int control;
			do{
				control=0;
				r= int(Nc*rnd.Rannyu());
				for(int k=0; k<j; k++){
					if(population[i][k]==r){control=1;}
				}
			}
			while(control==1);
			population[i][j]=r;
		}
	}
	rnd.SaveSeed();
}

void GA :: Distanze(){						//calcolo distanze salvate in dist
	dist.resize(Np);
	for (int i=0; i<Np; i++){
		dist[i]=0;
		for (int j=0; j<Nc; j++){
			int c1=population[i][j];
			int c2=population[i][(j+1)%Nc];
			double d=sqrt(pow(cities[c1][0]-cities[c2][0],2)+pow(cities[c1][1]-cities[c2][1],2));
			dist[i]=dist[i]+d;
		}
	}
}

void GA :: Order(){						//gerarchia distanze in idist 
	idist.resize(Np);
	for (int i=0; i<Np; i++){
		int cont=0;
		for (int j=0; j<Np; j++){
			if (dist[i]>dist[j]){cont+=1;}
		}
		for (int k=0; k<i; k++){
			if (idist[k]==cont){cont+=1;}
		}
		idist[i]=cont;
	}	
}

int GA :: Selection(){						//selezione dei migliori
	int r=int(Np*pow(rnd.Rannyu(),2));
	return r;
}

void GA:: Crossover() {
	int i=0;
	while(i<Np){
		int im=Selection();
		int ip=Selection();
		int imama=0.;
		int ipapa=0.;
		for (int k=0; k<Np; k++){
			if (idist[k]==im){imama=k;}
			if (idist[k]==ip){ipapa=k;}
		}
		if (rnd.Rannyu()<0.5){				//mama e papa copiati uguali
			for(int j=0; j<Nc; j++){
				newpopulation[i][j]=population[imama][j];
				newpopulation[i+1][j]=population[ipapa][j];
			}
			i=i+2;
		}
		else{						//crossover
			int cut=int(Nc*rnd.Rannyu());
			for(int j=0; j<cut; j++){
				newpopulation[i][j]=population[imama][j];
				newpopulation[i+1][j]=population[ipapa][j];
			}
			for(int j=cut; j<Nc; j++){
				int papa;
				int mama;
				int control;
				int cont=0;
				do {
					control=0;
					papa=population[ipapa][cont];
					cont=cont+1;
					for(int k=0; k<j; k++){
						if (papa==newpopulation[i][k]) {control=1;}
					}
				}
				while(control==1);
				newpopulation[i][j]=papa;
				cont=0;
				do {
					control=0;
					mama=population[imama][cont];
					cont=cont+1;
					for(int k=0; k<j; k++){
						if (mama==newpopulation[i+1][k]) {control=1;}
					}
				}
				while(control==1);
				newpopulation[i+1][j]=mama;
			}
			i=i+2;
		}
	}
}

void GA :: Mutation(int type=0){
	if (type==0){					//permutazione di due città
		for (int i=0; i<Np; i++){
			if (rnd.Rannyu()<0.1){
				int a=int((Nc-1)*rnd.Rannyu())+1;
				int b=int((Nc-1)*rnd.Rannyu())+1;
				int c=newpopulation[i][a];
				newpopulation[i][a]=newpopulation[i][b];
				newpopulation[i][b]=c;
			}
		}
	}
	else if (type==1){					//shift di ncit città di npos posizioni a partire dalla città ic
		for (int i=0; i<Np; i++){
			if (rnd.Rannyu()<0.05){
				vector<int> temp(Nc-1);
				for (int j=0; j<(Nc-1); j++){temp[j]=newpopulation[i][j+1];}
				int ncit=(Nc-1)*rnd.Rannyu();	
				int npos=(Nc-1-ncit)*rnd.Rannyu();
				int ic=(Nc-1)*rnd.Rannyu();
				vector<int> t_cit(ncit);
				vector<int> t_pos(npos);
				for (int j=0; j<ncit; j++){t_cit[j]=temp[(ic+j)%(Nc-1)];}
				for (int j=0; j<npos; j++){t_pos[j]=temp[(ic+ncit+j)%(Nc-1)];}
				for (int j=0; j<npos; j++){temp[(ic+j)%(Nc-1)]=t_pos[j];}
				for (int j=0; j<ncit; j++){temp[(ic+npos+j)%(Nc-1)]=t_cit[j];}
				for (int j=0; j<(Nc-1); j++){newpopulation[i][j+1]=temp[j];}
			}
		}
	}
	if (type==2){						//scambio di ncit adiacenti
		for (int i=0; i<Np; i++){
			if (rnd.Rannyu()<0.05){
				vector<int> temp(Nc-1);
				for (int j=0; j<(Nc-1); j++){temp[j]=newpopulation[i][j+1];}
				int ncit=((Nc-1)/2)*rnd.Rannyu();	
				int ic=(Nc-1)*rnd.Rannyu();
				vector<int> t_cit(ncit);
				for (int j=0; j<ncit; j++){t_cit[j]=temp[(ic+j)%(Nc-1)];}
				for (int j=0; j<ncit; j++){temp[(ic+j)%(Nc-1)]=temp[(ic+ncit+j)%(Nc-1)];}
				for (int j=0; j<ncit; j++){temp[(ic+ncit+j)%(Nc-1)]=t_cit[j];}
				for (int j=0; j<(Nc-1); j++){newpopulation[i][j+1]=temp[j];}
			}
		}
	}
	if (type==3){						//inversione ordine di ncit adiacenti
		for (int i=0; i<Np; i++){
			if (rnd.Rannyu()<0.1){
				vector<int> temp(Nc-1);
				for (int j=0; j<(Nc-1); j++){temp[j]=newpopulation[i][j+1];}
				int ncit=(Nc-1)*rnd.Rannyu();	
				int ic=(Nc-1)*rnd.Rannyu();
				vector<int> t_cit(ncit);
				for (int j=0; j<ncit; j++){t_cit[j]=temp[(ic+j)%(Nc-1)];}
				for (int j=0; j<ncit; j++){temp[(ic+j)%(Nc-1)]=t_cit[(ncit-1-j)];}
				for (int j=0; j<(Nc-1); j++){newpopulation[i][j+1]=temp[j];}
			}
		}
	}

}

void GA :: NewGeneration(){
	Distanze();
	Order();
	Crossover();
	Mutation(0);
	Mutation(1);
	Mutation(2);
	Mutation(3);
	for (int i=0; i<Np; i++){
		for (int j=0; j<Nc; j++){
			population[i][j]=newpopulation[i][j];
		}
	}
} 

void GA :: NLBest(int n){
	nL_best.resize(n);
	Distanze();
	Order();
	int ibest;
	for (int i=0; i<n; i++){
		for (int j=0; j<Np; j++){
			if (i==idist[j]){ibest=j;}
		}
		nL_best[i]=dist[ibest];
	}
}

int GA :: iBest(){
	Distanze();
	Order();
	int ibest;
	for (int i=0; i<Np; i++){
		if(0==idist[i]){ibest=i;}
	}
	return ibest;
}
















