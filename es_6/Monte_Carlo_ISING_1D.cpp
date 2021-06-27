/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

void stampa_file (const char *nome,double v[], int N){
ofstream x;
x.open(nome);
if (x.is_open()){
	for (int k=0; k<N; k++){x<<v[k]<<endl;}	
	}else cerr<<"problema: non si è aperto il file "<<nome<<endl;
x.close();
}

//media delle prime N componenti di un vettore v
double media (double v[], int N) {
	double sum=0;
	for (int i=0; i<N; i++){
		sum=sum+v[i];
	}
	return sum/N;
}

//varianza delle prime N componenti di un vettore v
double varianza (double v[], int N){
	double v1=0;
	double v2=0;
	for (int i=0; i<N; i++){
		v1=v1+v[i];
		v2=v2+v[i]*v[i];
	}
	return (v2/N)-((v1*v1)/(N*N));
}

int main(){
remove("output.ene.0");
remove("output.heat.0"); 
remove("output.mag.0"); 
remove("output.chi.0"); 

/*
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
*/

// EQUILIBRIO E CORRELAZIONE 

int n_blocchi=1000; 
int n_componenti=1;
double energy_eq[n_blocchi];

Input("input.dat"); //Inizialization
for (int i=0; i<n_blocchi; i++){
	Reset(i);
	for(int j=0; j<n_componenti; j++){
		Move(1);
		Measure();
	}
      	Accumulate();
      	Averages(i);
      	energy_eq[i]=stima_u;
}
stampa_file("energy_eq_met",energy_eq,n_blocchi);

Input("input.dat"); //Inizialization
for (int i=0; i<n_blocchi; i++){
	Reset(i);
	for(int j=0; j<n_componenti; j++){
		Move(0);
		Measure();
	}
      	Accumulate();
      	Averages(i);
      	energy_eq[i]=stima_u;
}
stampa_file("energy_eq_gib",energy_eq,n_blocchi);

//SIMULAZIONE
n_blocchi=500; 
n_componenti=200;
double en[n_blocchi];
double ct[n_blocchi];
double ma[n_blocchi];
double sm[n_blocchi];
double en_med[16];
double ct_med[16];
double ma_med[16];
double sm_med[16];
double en_err[16];
double ct_err[16];
double ma_err[16];
double sm_err[16];
Input("input.dat"); //Inizialization

//METROPOLIS metro=1
for (int k=0; k<16; k++){
	temp=0.5+k*0.1;
	for (int i=0; i<200; i++){Move(1);}  //equilibrazione
	for (int i=0; i<n_blocchi; i++){
		Reset(i);
		for(int j=1; j<=n_componenti; j++){
			Move(1);
			Measure();
		}
      		Accumulate();
      		Averages(i);
      		en[i]=stima_u;
      		ct[i]=stima_c;
      		ma[i]=stima_m;
      		sm[i]=stima_x;
	}
	en_med[k]=media(en,n_blocchi);
	ct_med[k]=media(ct,n_blocchi);
	sm_med[k]=media(sm,n_blocchi);
	en_err[k]=pow(varianza(en,n_blocchi)/n_blocchi,0.5);
	ct_err[k]=pow(varianza(ct,n_blocchi)/n_blocchi,0.5);
	sm_err[k]=pow(varianza(sm,n_blocchi)/n_blocchi,0.5);
}

stampa_file("en_M", en_med, 16);
stampa_file("ct_M", ct_med, 16);
stampa_file("sm_M", sm_med, 16);
stampa_file("en_err_M", en_err, 16);
stampa_file("ct_err_M", ct_err, 16);
stampa_file("sm_err_M", sm_err, 16);

//METROPOLIS metro=1
for (int k=0; k<16; k++){
	temp=0.5+k*0.1;
	for (int i=0; i<200; i++){Move(0);}  //equilibrazione
	for (int i=0; i<n_blocchi; i++){
		Reset(i);
		for(int j=1; j<=n_componenti; j++){
			Move(0);
			Measure();
		}
      		Accumulate();
      		Averages(i);
      		en[i]=stima_u;
      		ct[i]=stima_c;
      		ma[i]=stima_m;
      		sm[i]=stima_x;
	}
	en_med[k]=media(en,n_blocchi);
	ct_med[k]=media(ct,n_blocchi);
	sm_med[k]=media(sm,n_blocchi);
	en_err[k]=pow(varianza(en,n_blocchi)/n_blocchi,0.5);
	ct_err[k]=pow(varianza(ct,n_blocchi)/n_blocchi,0.5);
	sm_err[k]=pow(varianza(sm,n_blocchi)/n_blocchi,0.5);
}

stampa_file("en_G", en_med, 16);
stampa_file("ct_G", ct_med, 16);
stampa_file("sm_G", sm_med, 16);
stampa_file("en_err_G", en_err, 16);
stampa_file("ct_err_G", ct_err, 16);
stampa_file("sm_err_G", sm_err, 16);

h=0.02;

for (int k=0; k<16; k++){
	temp=0.5+k*0.1;
	for (int i=0; i<200; i++){Move(1);}  //equilibrazione
	for (int i=0; i<n_blocchi; i++){
		Reset(i);
		for(int j=1; j<=n_componenti; j++){
			Move(1);
			Measure();
		}
      		Accumulate();
      		Averages(i);
      		en[i]=stima_u;
      		ct[i]=stima_c;
      		ma[i]=stima_m;
      		sm[i]=stima_x;
	}
	ma_med[k]=media(ma,n_blocchi);
	ma_err[k]=pow(varianza(ma,n_blocchi)/n_blocchi,0.5);
}
stampa_file("ma_M", ma_med, 16);
stampa_file("ma_err_M", ma_err, 16);

for (int k=0; k<16; k++){
	temp=0.5+k*0.1;
	for (int i=0; i<200; i++){Move(0);}  //equilibrazione
	for (int i=0; i<n_blocchi; i++){
		Reset(i);
		for(int j=1; j<=n_componenti; j++){
			Move(0);
			Measure();
		}
      		Accumulate();
      		Averages(i);
      		en[i]=stima_u;
      		ct[i]=stima_c;
      		ma[i]=stima_m;
      		sm[i]=stima_x;
	}
	ma_med[k]=media(ma,n_blocchi);
	ma_err[k]=pow(varianza(ma,n_blocchi)/n_blocchi,0.5);
}
stampa_file("ma_G", ma_med, 16);
stampa_file("ma_err_G", ma_err, 16);


  return 0;
}


void Input(const char *nome)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open(nome);

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}

void Input2(const char *nome)
{
  ifstream ReadInput;
  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;
//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();
   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();  
//Read input informations
  ReadInput.open(nome);
  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;
  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;
  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;
  ReadInput >> h;
  cout << "External field = " << h << endl << endl;   
  ReadInput >> metro; // if=1 Metropolis else Gibbs
  ReadInput >> nblk;
  ReadInput >> nstep;
  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();
//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
  n_props = 4; //Number of observables
//initial configuration
  ifstream read;
  read.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    read>>s[i];
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}




void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down, p_up;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    
    if(metro==1){ //Metropolis
    	attempted++;
	sm=s[o];
	energy_old=Boltzmann(sm, o);
	sm=-1*sm;
	energy_new=Boltzmann(sm,o);
	p=exp(-(energy_new-energy_old)/temp);
	if (p>1){p=1;}
	if (rnd.Rannyu()<p){
		s[o]=sm;
		accepted++;
	}
    }
    else{ //Gibbs sampling
	energy_up=exp(-Boltzmann(1,o)/temp);
	energy_down=exp(-Boltzmann(-1,o)/temp);
	p_up=energy_up/(energy_up+energy_down);
	if (rnd.Rannyu()<p_up){s[o]=+1;}
	else {s[o]=-1;}
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, c=0.0, m = 0.0, x=0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     c += pow(-J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]),2);   //<H^2> non proprio c
     m +=s[i];
  }
  x = pow (m,2);
  walker[iu] = u;
  walker[ic] = c;
  walker[im] = m;
  walker[ix] = x;
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Ene.open("output.ene.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();
    
    Heat.open("output.heat.0",ios::app);
    stima_c = pow(temp,-2)*((blk_av[ic]/blk_norm/(double)nspin)-pow(stima_u,2));
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();
    
    Mag.open("output.mag.0",ios::app);
    stima_m = blk_av[im]/blk_norm/(double)nspin;
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Mag.close();
    
    Chi.open("output.chi.0",ios::app);
    stima_x = pow(temp,-1)*(blk_av[ix]/blk_norm/(double)nspin);
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Chi.close();
    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
