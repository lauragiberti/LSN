/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <cstdio>
#include "MolDyn_NVE.h"

using namespace std;
void stampa(double v[], int N){
	for (int i=0; i<N; i++){
		cout<<v[i]<<endl;
	}
}

//varianza delle prime N componenti di un vettore v
double varianza (double v[], int N){
	double v1=0;
	double v2=0;
	for (int i=0; i<N; i++){
		v1=v1+v[i];
		v2=v2+v[i]*v[i];
	}
	return ((v2/N)-((v1*v1)/(N*N)))/N;
}

//media delle prime N componenti di un vettore v
double media (double v[], int N) {
	double sum=0;
	for (int i=0; i<N; i++){
		sum=sum+v[i];
	}
	return sum/N;
}

//salva in "risultati" la media a blocchi
//N=elementi di v;  n=componenti di ogni blocco;  N/n=numero di blocchi
void media_blocchi (double v[], double risultati[], int N, int n) {
	if (N%n!=0){cout<<"errore Nblocchi non divisibili"<<endl;}
	double medie_blocchi[N/n];
	int k=0;
	for (int i=0; i<(N/n); i++){
		double sum=0;
		for (int j=0; j<n; j++){
			sum=sum+v[k];
			k++;
		}
		medie_blocchi[i]=sum/n;
		risultati[i]=medie_blocchi[i];
	}
}

void stampa_file (const char *nome,double v[], int N){
ofstream x;
x.open(nome);
if (x.is_open()){
	for (int k=0; k<N; k++){x<<v[k]<<endl;}	
	}else cerr<<"problema: non si è aperto il file "<<nome<<endl;
x.close();
}


int main(){ 
remove("output_ekin.dat");
remove("output_epot.dat");
remove("output_etot.dat");
remove("output_temp.dat");

int n=10000;
double T[n/10];
double Etot[n/10];
double Ekin[n/10];
double Epot[n/10];

//esercizio1
Input();  
for(int i=1; i <= n; i++){
	Move();           //Move particles with Verlet algorithm
	if(i%iprint == 0) cout << "Numero passi: " << i << " di " << n << endl;
	if(i%10 == 0){
		Measure();     //Properties measurement
		T[(i/10)-1]=stima_temp;
		Etot[(i/10)-1]=stima_etot;
		Ekin[(i/10)-1]=stima_kin;
		Epot[(i/10)-1]=stima_pot;
	}
}
ConfFinal();         //Write final configuration to restart
stampa_file("T_init",T,n/10);
stampa_file("Etot_init",Etot,n/10);
stampa_file("Ekin_init",Ekin,n/10);
stampa_file("Epot_init",Epot,n/10);

Input2();
for (int i=1; i<=n; i++){
	Move();
	if(i%iprint == 0) cout << "Numero passi: " << i << " di " << n << endl;
	if(i%10==0){
		Measure();
		T[(i/10)-1]=stima_temp;
		Etot[(i/10)-1]=stima_etot;
		Ekin[(i/10)-1]=stima_kin;
		Epot[(i/10)-1]=stima_pot;
	}
}
stampa_file("T_nonrescale",T,n/10);
stampa_file("Etot_nonrescale",Etot,n/10);
stampa_file("Ekin_nonrescale",Ekin,n/10);
stampa_file("Epot_nonrescale",Epot,n/10);

Input2();
for (int i=0; i<5000; i++){
	Move();
	if (i%10 == 0){Rescale_v();}
}
ConfFinal();         //Write final configuration to restart

Input2();
for (int i=1; i<=n; i++){
	Move();
	if(i%iprint == 0) cout << "Numero passi: " << i << " di " << n << endl;
	if (i%10==0){
		Measure();
		T[(i/10)-1]=stima_temp;
		Etot[(i/10)-1]=stima_etot;
		Ekin[(i/10)-1]=stima_kin;
		Epot[(i/10)-1]=stima_pot;
	}
}
stampa_file("T_rescale",T,n/10);
stampa_file("Etot_rescale",Etot,n/10);
stampa_file("Ekin_rescale",Ekin,n/10);
stampa_file("Epot_rescale",Epot,n/10);

for (int i=0; i<(n/10); i++){
	Etot[i]=Etot[i]/npart;
	Ekin[i]=Ekin[i]/npart;
	Epot[i]=Epot[i]/npart;
}

//esercizio 2
double T_ave[n/100];
double Etot_ave[n/100];
double Ekin_ave[n/100];
double Epot_ave[n/100];


media_blocchi(T, T_ave, n/10, 10);
media_blocchi(Etot, Etot_ave, n/10, 10);
media_blocchi(Ekin, Ekin_ave, n/10, 10);
media_blocchi(Epot, Epot_ave, n/10, 10);


double T_ave_smooth[n/100];
double Etot_ave_smooth[n/100];
double Ekin_ave_smooth[n/100];
double Epot_ave_smooth[n/100];

double T_ave_smooth_var[n/100];
double Etot_ave_smooth_var[n/100];
double Ekin_ave_smooth_var[n/100];
double Epot_ave_smooth_var[n/100];

for (int i=0; i<n/100; i++){
	T_ave_smooth[i]=media(T_ave, i+1);
	Etot_ave_smooth[i]=media(Etot_ave, i+1);
	Ekin_ave_smooth[i]=media(Ekin_ave, i+1);
	Epot_ave_smooth[i]=media(Epot_ave, i+1);
	
	T_ave_smooth_var[i]=varianza(T_ave, i+1);
	Etot_ave_smooth_var[i]=varianza(Etot_ave, i+1);
	Ekin_ave_smooth_var[i]=varianza(Ekin_ave, i+1);
	Epot_ave_smooth_var[i]=varianza(Epot_ave, i+1);
}

stampa_file("T_smooth",T_ave_smooth,n/100);
stampa_file("Etot_smooth",Etot_ave_smooth,n/100);
stampa_file("Ekin_smooth",Ekin_ave_smooth,n/100);
stampa_file("Epot_smooth",Epot_ave_smooth,n/100);

stampa_file("T_smooth_var",T_ave_smooth_var,n/100);
stampa_file("Etot_smooth_var",Etot_ave_smooth_var,n/100);
stampa_file("Ekin_smooth_var",Ekin_ave_smooth_var,n/100);
stampa_file("Epot_smooth_var",Epot_ave_smooth_var,n/100);


return 0;
}



void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
   return;
}

void Input2(void){ //prepara tutto per una ripartenza
  ifstream ReadInput,ReadConf, ReadOld;
  double ep, ek, pr, et, vir;
  
  ReadInput.open("input.dat"); //Read input
  ReadInput >> temp;
  ReadInput >> npart;
  ReadInput >> rho;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Read initial configuration
  cout << "Read initial configuration from file config.final e old.final" << endl << endl;
  ReadConf.open("config.final");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();
  ReadOld.open("old.final");
  for (int i=0; i<npart; ++i){
    ReadOld >> xold[i] >> yold[i] >> zold[i];
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
  }
  ReadOld.close();

   return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "configurazione finale e precedente salvate in config.final e old.final" << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  
  WriteConf.open("old.final");
  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}


void Rescale_v(void){
  Move();
  for (int i=0; i<npart;i++){
    vx[i] = Pbc(x[i] - xold[i])/(delta);
    vy[i] = Pbc(y[i] - yold[i])/(delta);
    vz[i] = Pbc(z[i] - zold[i])/(delta);
  }
  
  double v2_mis=0;
  double fs;
  for (int i=0; i<npart; i++){
    v2_mis=v2_mis+vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];
  }
  v2_mis=v2_mis/npart;
  fs=sqrt(3*temp/v2_mis);
  
  for (int i=0; i<npart;i++){
    vx[i] = vx[i]*fs;
    vy[i] = vy[i]*fs;
    vz[i] = vz[i]*fs;
  }
  
 for (int i=0; i<npart;i++){
    xold[i] = Pbc(-vx[i]*delta + x[i]);
    yold[i] = Pbc(-vy[i]*delta + y[i]);
    zold[i] = Pbc(-vz[i]*delta + z[i]);
  }
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
