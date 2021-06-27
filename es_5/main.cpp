#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

void stampa(double v[], int N){
	for (int i=0; i<N; i++){
		cout<<v[i]<<endl;
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

double media (double v[], int N) {
	double sum=0;
	for (int i=0; i<N; i++){
		sum=sum+v[i];
	}
	return sum/N;
}

double varianza (double v[], int N){
	double v1=0;
	double v2=0;
	for (int i=0; i<N; i++){
		v1=v1+v[i];
		v2=v2+v[i]*v[i];
	}
	return (v2/N)-((v1*v1)/(N*N));
}

double phi1 (double x,double y,double z){
	double phi=pow(2/M_PI,3/2)*exp(-2*sqrt(x*x+y*y+z*z));
	return phi;
}

double phi2 (double x,double y,double z){
	double phi=z*z*exp(-sqrt(x*x+y*y+z*z));
	return phi;
}
 
int main (int argc, char *argv[]){

//NUMERI CASUALI
   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
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
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

double l=1.2;
double x=0., y=0., z=1.;
double xnew, ynew, znew;
double p;
int step=1000000;
int blocchi=1000;
double r;
double sum=0;
double r_ave[blocchi];
double r_ave_smooth[blocchi];
double r_var_smooth[blocchi];
ofstream prx;
ofstream pry;
ofstream prz;





//UNIFORME PHI1
int k=0;
prx.open("x_unif_1");
pry.open("y_unif_1");
prz.open("z_unif_1");

for (int i=0; i<step/10; i++){
	xnew=x+(rnd.Rannyu()*2*l-l);
		ynew=y+(rnd.Rannyu()*2*l-l);
		znew=z+(rnd.Rannyu()*2*l-l);
		p=phi1(xnew,ynew,znew)/phi1(x,y,z);
		if (p>1){p=1;}
		double rand=rnd.Rannyu();
		if(rand<p){
			x=xnew;
			y=ynew;
			z=znew;
		}
}

for (int i=0; i<blocchi; i++){
	sum=0;
	for(int j=0; j<(step/blocchi); j++){
		xnew=x+(rnd.Rannyu()*2*l-l);
		ynew=y+(rnd.Rannyu()*2*l-l);
		znew=z+(rnd.Rannyu()*2*l-l);
		p=phi1(xnew,ynew,znew)/phi1(x,y,z);
		if (p>1){p=1;}
		double rand=rnd.Rannyu();
		if(rand<p){
			x=xnew;
			y=ynew;
			z=znew;
			k++;
		}
		r=sqrt(x*x+y*y+z*z);
		sum=sum+r;
		if (j%100==0){
			prx<<x<<endl;
			pry<<y<<endl;
			prz<<z<<endl;
		}
	}
	r_ave[i]=sum/(step/blocchi);
}
prx.close();
pry.close();
prz.close();

cout<<"UNIFORME PHI 1"<<endl;
cout<<"probabilità accettazzione: "<<k*100/step<<"%"<<endl;
cout<<"media: "<<media(r_ave, step/blocchi)<<endl;
cout<<"err: "<<sqrt(varianza(r_ave, step/blocchi)/(step/blocchi))<<endl<<endl;

for (int i=0; i<blocchi; i++){
	r_ave_smooth[i]=media(r_ave, i+1);
	r_var_smooth[i]=varianza(r_ave, i+1)/(i+1);
}

stampa_file("r_unif_phi1",r_ave_smooth, blocchi);
stampa_file("var_unif_phi1",r_var_smooth, blocchi);

//UNIFORME PHI2
k=0;
l=2.9;
x=0., y=0.; z=1.;
prx.open("x_unif_2");
pry.open("y_unif_2");
prz.open("z_unif_2");

for (int i=0; i<step/10; i++){
	xnew=x+(rnd.Rannyu()*2*l-l);
		ynew=y+(rnd.Rannyu()*2*l-l);
		znew=z+(rnd.Rannyu()*2*l-l);
		p=phi1(xnew,ynew,znew)/phi1(x,y,z);
		if (p>1){p=1;}
		double rand=rnd.Rannyu();
		if(rand<p){
			x=xnew;
			y=ynew;
			z=znew;
		}
}

for (int i=0; i<blocchi; i++){
	sum=0;
	for(int j=0; j<(step/blocchi); j++){
		xnew=x+(rnd.Rannyu()*2*l-l);
		ynew=y+(rnd.Rannyu()*2*l-l);
		znew=z+(rnd.Rannyu()*2*l-l);
		p=phi2(xnew,ynew,znew)/phi2(x,y,z);
		if (p>1){p=1;}
		double rand=rnd.Rannyu();
		if(rand<p){
			x=xnew;
			y=ynew;
			z=znew;
			k++;
		}
		r=sqrt(x*x+y*y+z*z);
		sum=sum+r;
		if (j%100==0){
			prx<<x<<endl;
			pry<<y<<endl;
			prz<<z<<endl;
		}
	}
	r_ave[i]=sum/(step/blocchi);
}
prx.close();
pry.close();
prz.close();
cout<<"UNIFORME PHI 2"<<endl;
cout<<"probabilità accettazzione: "<<k*100/step<<"%"<<endl;
cout<<"media: "<<media(r_ave, step/blocchi)<<endl;
cout<<"err: "<<sqrt(varianza(r_ave, step/blocchi)/(step/blocchi))<<endl<<endl;

for (int i=0; i<blocchi; i++){
	r_ave_smooth[i]=media(r_ave, i+1);
	r_var_smooth[i]=varianza(r_ave, i+1)/(i+1);
}

stampa_file("r_unif_phi2",r_ave_smooth, blocchi);
stampa_file("var_unif_phi2",r_var_smooth, blocchi);


//GAUSSIANA PHI1

k=0;
l=0.75;
x=0., y=0.; z=1.;

prx.open("x_gauss_1");
pry.open("y_gauss_1");
prz.open("z_gauss_1");

for (int i=0; i<step/10; i++){
	xnew=x+(rnd.Rannyu()*2*l-l);
		ynew=y+(rnd.Rannyu()*2*l-l);
		znew=z+(rnd.Rannyu()*2*l-l);
		p=phi1(xnew,ynew,znew)/phi1(x,y,z);
		if (p>1){p=1;}
		double rand=rnd.Rannyu();
		if(rand<p){
			x=xnew;
			y=ynew;
			z=znew;
		}
}

for (int i=0; i<blocchi; i++){
	sum=0;
	for(int j=0; j<(step/blocchi); j++){
		xnew=rnd.Gauss(x,l);
		ynew=rnd.Gauss(y,l);
		znew=rnd.Gauss(z,l);
		p=phi1(xnew,ynew,znew)/phi1(x,y,z);
		if (p>1){p=1;}
		double rand=rnd.Rannyu();
		if(rand<p){
			x=xnew;
			y=ynew;
			z=znew;
			k++;
		}
		r=sqrt(x*x+y*y+z*z);
		sum=sum+r;
		if (j%100==0){
			prx<<x<<endl;
			pry<<y<<endl;
			prz<<z<<endl;
		}
	}
	r_ave[i]=sum/(step/blocchi);
}
prx.close();
pry.close();
prz.close();
cout<<"GAUSS PHI 1"<<endl;
cout<<"probabilità accettazzione: "<<k*100/step<<"%"<<endl;
cout<<"media: "<<media(r_ave, step/blocchi)<<endl;
cout<<"err: "<<sqrt(varianza(r_ave, step/blocchi)/(step/blocchi))<<endl<<endl;

for (int i=0; i<blocchi; i++){
	r_ave_smooth[i]=media(r_ave, i+1);
	r_var_smooth[i]=varianza(r_ave, i+1)/(i+1);
}

stampa_file("r_gauss_phi1",r_ave_smooth, blocchi);
stampa_file("var_gauss_phi1",r_var_smooth, blocchi);

//GAUSS PHI2
k=0;
l=1.85;
x=0., y=0.; z=1.;
prx.open("x_gauss_2");
pry.open("y_gauss_2");
prz.open("z_gauss_2");

for (int i=0; i<step/10; i++){
	xnew=x+(rnd.Rannyu()*2*l-l);
		ynew=y+(rnd.Rannyu()*2*l-l);
		znew=z+(rnd.Rannyu()*2*l-l);
		p=phi1(xnew,ynew,znew)/phi1(x,y,z);
		if (p>1){p=1;}
		double rand=rnd.Rannyu();
		if(rand<p){
			x=xnew;
			y=ynew;
			z=znew;
		}
}

for (int i=0; i<blocchi; i++){
	sum=0;
	for(int j=0; j<(step/blocchi); j++){
		xnew=rnd.Gauss(x,l);
		ynew=rnd.Gauss(y,l);
		znew=rnd.Gauss(z,l);
		p=phi2(xnew,ynew,znew)/phi2(x,y,z);
		if (p>1){p=1;}
		double rand=rnd.Rannyu();
		if(rand<p){
			x=xnew;
			y=ynew;
			z=znew;
			k++;
		}
		r=sqrt(x*x+y*y+z*z);
		sum=sum+r;
		if (j%100==0){
			prx<<x<<endl;
			pry<<y<<endl;
			prz<<z<<endl;
		}
	}
	r_ave[i]=sum/(step/blocchi);
}
prx.close();
pry.close();
prz.close();

cout<<"UNIFORME PHI 2"<<endl;
cout<<"probabilità accettazzione: "<<k*100/step<<"%"<<endl;
cout<<"media: "<<media(r_ave, step/blocchi)<<endl;
cout<<"err: "<<sqrt(varianza(r_ave, step/blocchi)/(step/blocchi))<<endl<<endl;

for (int i=0; i<blocchi; i++){
	r_ave_smooth[i]=media(r_ave, i+1);
	r_var_smooth[i]=varianza(r_ave, i+1)/(i+1);
}

stampa_file("r_gauss_phi2",r_ave_smooth, blocchi);
stampa_file("var_gauss_phi2",r_var_smooth, blocchi);

   rnd.SaveSeed();
   return 0;

}

