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

double varianza (double v[], int N){
	double v1=0;
	double v2=0;
	for (int i=0; i<N; i++){
		v1=v1+v[i];
		v2=v2+v[i]*v[i];
	}
	return ((v2/N)-((v1*v1)/(N*N)))/N;
}

double media (double v[], int N) {
	double sum=0;
	for (int i=0; i<N; i++){
		sum=sum+v[i];
	}
	return sum/N;
}

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

int n_blocchi=100;
int n_componenti=200;

double gr_ave[100]={0}, gr_ave2[100]={0}, gr_err[100];
double etot[n_blocchi], pres[n_blocchi];
double etot_smooth[n_blocchi], pres_smooth[n_blocchi];
double etot_err_smooth[n_blocchi], pres_err_smooth[n_blocchi];


//ARGON SOLIDO
Input("input.solid");
for(int i=0; i < 1000; i++){
	Move();
	if (i%10 == 0){Rescale_v();}
}
ConfFinal();

Input2("input.solid");
for(int i=0; i<n_blocchi; i++){
	Reset();
	cout << i << " di " << n_blocchi << endl;
	for(int j=0; j<n_componenti; j++){
		Move();
		Measure();
	}
	Averages();
	etot[i]=stima_etot;
	pres[i]=stima_pres;
	for (int i=0; i<100; i++){
		gr_ave[i]+= gr[i];
		gr_ave2[i]+= gr[i]*gr[i];
	}	
}

for (int i=0; i<n_blocchi; i++){
	etot_smooth[i]=media(etot, i+1);
	etot_err_smooth[i]=pow(varianza(etot, i+1)/n_blocchi,0.5);
	pres_smooth[i]=media(pres, i+1);
	pres_err_smooth[i]=pow(varianza(pres, i+1)/n_blocchi,0.5);
}

for (int i=0; i<100; i++){
	gr_ave[i] = gr_ave[i]/n_blocchi;
	gr_ave2[i]= gr_ave2[i]/n_blocchi;
	gr_err[i]=pow((gr_ave2[i]-gr_ave[i]*gr_ave[i])/n_blocchi,0.5);
}

stampa_file("gr_S",gr_ave, 100);
stampa_file("gr_err_S",gr_err, 100);
stampa_file("etot_S", etot_smooth, n_blocchi);
stampa_file("etot_err_S", etot_err_smooth, n_blocchi);
stampa_file("pres_S", pres_smooth, n_blocchi);
stampa_file("pres_err_S", pres_err_smooth, n_blocchi);

//ARGON LIQUIDO
for (int i=0; i<100; i++){
	gr_ave[i]=0;
	gr_ave2[i]=0;
}

Input("input.liquid");
for(int i=0; i < 1000; i++){
	Move();
	if (i%10 == 0){Rescale_v();}
}
ConfFinal();

Input2("input.liquid");
for(int i=0; i<n_blocchi; i++){
	Reset();
	cout << i << " di " << n_blocchi << endl;
	for(int j=0; j<n_componenti; j++){
		Move();
		Measure();
	}
	Averages();
	etot[i]=stima_etot;
	pres[i]=stima_pres;
	for (int i=0; i<100; i++){
		gr_ave[i]+= gr[i];
		gr_ave2[i]+= gr[i]*gr[i];
	}	
}

for (int i=0; i<n_blocchi; i++){
	etot_smooth[i]=media(etot, i+1);
	etot_err_smooth[i]=pow(varianza(etot, i+1)/n_blocchi,0.5);
	pres_smooth[i]=media(pres, i+1);
	pres_err_smooth[i]=pow(varianza(pres, i+1)/n_blocchi,0.5);
}

for (int i=0; i<100; i++){
	gr_ave[i] = gr_ave[i]/n_blocchi;
	gr_ave2[i]= gr_ave2[i]/n_blocchi;
	gr_err[i]=pow((gr_ave2[i]-gr_ave[i]*gr_ave[i])/n_blocchi,0.5);
}

stampa_file("gr_L",gr_ave, 100);
stampa_file("gr_err_L",gr_err, 100);
stampa_file("etot_L", etot_smooth, n_blocchi);
stampa_file("etot_err_L", etot_err_smooth, n_blocchi);
stampa_file("pres_L", pres_smooth, n_blocchi);
stampa_file("pres_err_L", pres_err_smooth, n_blocchi);

//ARGON GAS
for (int i=0; i<100; i++){
	gr_ave[i]=0;
	gr_ave2[i]=0;
}

Input("input.gas");
for(int i=0; i < 1000; i++){
	Move();
	if (i%10 == 0){Rescale_v();}
}
ConfFinal();

Input2("input.gas");
for(int i=0; i<n_blocchi; i++){
	Reset();
	cout << i << " di " << n_blocchi << endl;
	for(int j=0; j<n_componenti; j++){
		Move();
		Measure();
	}
	Averages();
	etot[i]=stima_etot;
	pres[i]=stima_pres;
	for (int i=0; i<100; i++){
		gr_ave[i]+= gr[i];
		gr_ave2[i]+= gr[i]*gr[i];
	}	
}

for (int i=0; i<n_blocchi; i++){
	etot_smooth[i]=media(etot, i+1);
	etot_err_smooth[i]=pow(varianza(etot, i+1)/n_blocchi,0.5);
	pres_smooth[i]=media(pres, i+1);
	pres_err_smooth[i]=pow(varianza(pres, i+1)/n_blocchi,0.5);
}

for (int i=0; i<100; i++){
	gr_ave[i] = gr_ave[i]/n_blocchi;
	gr_ave2[i]= gr_ave2[i]/n_blocchi;
	gr_err[i]=pow((gr_ave2[i]-gr_ave[i]*gr_ave[i])/n_blocchi,0.5);
}

stampa_file("gr_G",gr_ave, 100);
stampa_file("gr_err_G",gr_err, 100);
stampa_file("etot_G", etot_smooth, n_blocchi);
stampa_file("etot_err_G", etot_err_smooth, n_blocchi);
stampa_file("pres_G", pres_smooth, n_blocchi);
stampa_file("pres_err_G", pres_err_smooth, n_blocchi);

return 0;
}



void Input(const char *nome){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open(nome); //Read input

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

void Input2(const char *nome){ //prepara tutto per una ripartenza
  ifstream ReadInput,ReadConf, ReadOld;
  double ep, ek, pr, et, vir;
  
  ReadInput.open(nome); //Read input
  ReadInput >> temp;
  ReadInput >> npart;
  ReadInput >> rho;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  bin_size = box/200.0;
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

void Reset(void){
//reset g(r)
  cont=0;
  stima_kin=0;
  stima_pot=0;
  stima_temp=0;
  stima_etot=0;
  stima_pres=0;
  for (int k=0; k<100; k++){gr[k]=0;}
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
  double v, w, t, vij, wij;
  double dx, dy, dz, dr;

  v = 0.0; //reset observables
  w=0.0;
  t = 0.0;
  cont=cont+1;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);
     
     // g(r)
     for (int k=0; k<100; k++){
     	if (dr>k*bin_size && dr<(k+1)*bin_size){
     		gr[k]=gr[k]+2;
     	}
     }

     if(dr < rcut){
       vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
       wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);
       v += vij;
       w += wij;
     }
     }
    } 
  v = 4.0 * v;
  w = 48.0 * w / 3.0;         

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
  stima_pot = stima_pot + v; //Potential energy per particle
  stima_kin = stima_kin + t; //Kinetic energy per particle
  stima_temp = stima_temp + t; //Temperature
  stima_etot = stima_etot + (t+v); //Total energy per particle
  stima_pres = stima_pres+w; //Pressure
    return;
}

void Averages(void){
  ofstream Epot, Ekin, Etot, Temp;
  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);

  stima_pot=stima_pot/cont/(double)npart;
  stima_kin=stima_kin/cont/(double)npart;
  stima_temp=(2.0 / 3.0) * stima_temp/cont/(double)npart;
  stima_etot=stima_etot/cont/(double)npart;
  stima_pres = rho*temp+stima_pres/vol/cont/(double)npart;
  for (int k=0; k<100; k++){
  double deltav=(4/3)*M_PI*(pow((0.005+0.01*k+0.01)*box,3)-pow((0.005+0.01*k)*box,3));
  gr[k]=gr[k]/cont/rho/m_part/deltav;
  }
  
  Epot << stima_pot  << endl;
  Ekin << stima_kin  << endl;
  Temp << stima_temp << endl;
  Etot << stima_etot << endl;

  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
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
