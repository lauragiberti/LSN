#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

//FUNZIONI

double funzione (double x){
	return (M_PI*cos(M_PI*x/2.)/2.);
}

double distribuzione (double x){
	return((3./2.)*(1.-x*x));
}

void stepRWcubo(double r, double v[], int l=1){
	if (r<0 && r>1){cout<<"errore numero non compreso tra [0,1]"<<endl;}
	if (r>0 && r<(1./6.)){v[0]=v[0]-l;}
	if (r>(1./6.) && r<(2./6.)){v[0]=v[0]+l;}
	if (r>(2./6.) && r<(3./6.)){v[1]=v[1]-l;}
	if (r>(3./6.) && r<(4./6.)){v[1]=v[1]+l;}
	if (r>(4./6.) && r<(5./6.)){v[2]=v[2]-l;}
	if (r>(5./6.) && r<1){v[2]=v[2]-l;}
}

void stepRWcont(double theta, double phi, double v[], int l=1){
	if (theta<0 && theta>2*M_PI){cout<<"errore theta non compreso tra [0,2PI]"<<endl;}
	if (phi<0 && phi>M_PI){cout<<"errore phi non compreso tra [0,PI]"<<endl;}
	v[0]=v[0]+l*sin(phi)*cos(theta);
	v[1]=v[1]+l*sin(phi)*sin(theta);
	v[2]=v[2]+l*cos(phi);
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


//somma vettoriale
void somma (double a[], double b[], double risultato[], int N){
	for (int i=0; i<N; i++){
		risultato[i]=a[i]+b[i];
	}
}

//stampa su terminale
void stampa(double v[], int N){
	for (int i=0; i<N; i++){
		cout<<v[i]<<endl;
	}
}

//stampa su file
void stampa_file (const char *nome,double v[], int N){
ofstream x;
x.open(nome);
if (x.is_open()){
	for (int k=0; k<N; k++){x<<v[k]<<endl;}	
	}else cerr<<"problema: non si è aperto il file "<<nome<<endl;
x.close();
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


/*
ESERCIZIO 1 parte 1
calcolo integrale con distribuzione uniforme
*/

int M=10000;
double I[M]={0};
double err[M]={0};
int hit=10000;

for (int j=0; j<M; j++){
	double somma=0;
	double somma2=0;
	for (int i=0; i<hit; i++){
		double x=rnd.Rannyu();
		somma=somma+funzione(x);
		somma2=somma2+funzione(x)*funzione(x);
	}
I[j]=somma/hit;
err[j]=somma2/hit-I[j]*I[j];
}
double err_I=media(err,10000);
cout<<"errore dell'integrale calcolato con la distribuzione uniforme: "<<err_I<<endl;

double media_b[100]={0};
media_blocchi(I, media_b, 10000, 100);

double media_smooth[100]={0};
double var_smooth[100]={0};
for (int i=0; i<100; i++){
	media_smooth[i]=media(media_b,(i+1));
	var_smooth[i]=varianza(media_b,(i+1))/(i+1);
}

stampa_file("I_unif_medie", media_smooth, 100);
stampa_file("I_unif_var", var_smooth, 100);

/*
ESERCIZIO 1 parte 2
calcolo integrale con importance sampling
*/


double pmax=(3./2.);

for (int j=0; j<M; j++){
	double somma=0;
	double somma2=0;
	for (int i=0; i<hit; i++){
		double x;
		double T=0;
		while (T==0){
			x=rnd.Rannyu();
			double r=rnd.Rannyu();
			if (r<(distribuzione(x)/pmax)){T=1;}
		}
		somma=somma+funzione(x)/distribuzione(x);
		somma2=somma2+pow(funzione(x)/distribuzione(x),2);
	}
I[j]=somma/hit;
err[j]=somma2/hit-I[j]*I[j];
}
err_I=media(err,M);
cout<<"errore dell'integrale calcolato con l'importance sampling: "<<err_I<<endl;

media_blocchi(I, media_b, 10000, 100);

for (int i=0; i<100; i++){
	media_smooth[i]=media(media_b,(i+1));
	var_smooth[i]=varianza(media_b,(i+1))/(i+1);
}

stampa_file("I_imps_medie", media_smooth, 100);
stampa_file("I_imps_var", var_smooth, 100);


/*
RANDOM WALK
*/

int l=1;
int step=100;
double sumr2[step]={0};
double sumr[step]={0};
double r2[step]={0};
double r[step]={0};
double r2_medio[step]={0};
double r_medio[step]={0};
double media[step]={0};
double errRW[step]={0};


for (int j=0; j<10000; j++){
	double v[3]={0};
	for (int i=1; i<step; i++){
		double z=rnd.Rannyu();
		stepRWcubo(z,v,l);
		r2[i]=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
		r[i]=pow(r2[i],0.5);
	}
	somma(sumr2,r2,sumr2,step);
	somma(sumr,r,sumr,step);
}

for (int i=0; i<100;i++){
	r2_medio[i]=sumr2[i]/10000;
	r_medio[i]=sumr[i]/10000;
}

for (int i=0; i<100; i++){
	media[i]=r_medio[i];
	errRW[i]=sqrt(r2_medio[i]-r_medio[i]*r_medio[i]);
}

stampa_file("sigmaRW_cubo",errRW, 100);
stampa_file("mediaRW_cubo",media, 100);




for (int i=0;i<step; i++){
	sumr2[i]=0;
	sumr[i]=0;
}


for (int j=0; j<10000; j++){
	double v[3]={0};
	for (int i=1; i<step; i++){
		double theta=2*M_PI*rnd.Rannyu();
		double phi=M_PI*rnd.Rannyu();
		stepRWcont(theta,phi,v,l);
		r2[i]=v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
		r[i]=pow(r2[i],0.5);
	}
	somma(sumr2,r2,sumr2,step);
	somma(sumr,r,sumr,step);
	
}

for (int i=0; i<100;i++){
	r2_medio[i]=sumr2[i]/10000;
	r_medio[i]=sumr[i]/10000;
}

for (int i=0; i<100; i++){
	media[i]=r_medio[i];
	errRW[i]=sqrt(r2_medio[i]-r_medio[i]*r_medio[i]);
}
stampa_file("sigmaRW_cont",errRW, 100);
stampa_file("mediaRW_cont",media, 100);


return 0;

}
