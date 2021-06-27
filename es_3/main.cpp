#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"

using namespace std;

//FUNZIONI

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

//somma vettoriale
void somma (double a[], double b[], double risultato[], int N){
	for (int i=0; i<N; i++){
		risultato[i]=a[i]+b[i];
	}
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

//massimo
double max(double v[], int N){
	double M=v[0];
	for (int i=1; i<N; i++){
		if (v[i]>M){M=v[i];}
	}
	return M;
}

//minimo
double min(double v[], int N){
	double m=v[0];
	for (int i=1; i<N; i++){
		if (v[i]<m){m=v[i];}
	}
	return m;
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


double S0=100.;
double T=1.;
double K=100.;
double r=0.1;
double sigma=0.25;

double ST;
double C0[10000];
double P0[10000];
double x;

//SAMPLING DIRETTO

for (int i=0; i<10000; i++){
	ST=S0*exp((r-sigma*sigma/2)*T+sigma*rnd.Gauss(0,T));

	if ((ST-K)>0.){x=(ST-K);}
	else{x=0;}
	C0[i]=exp(-r*T)*x;

	if ((K-ST)>0.){x=(K-ST);}
	else{x=0;}
	P0[i]=exp(-r*T)*x;
}
cout<<"call: "<<media(C0,10000)<<"+-"<<pow(varianza(C0,10000),0.5)<<endl;
cout<<"put: "<<media(P0,10000)<<"+-"<<pow(varianza(P0,10000),0.5)<<endl;

double C0_blocchi[100];
double P0_blocchi[100];

double C0_smooth[100];
double C0_smooth_var[100];
double P0_smooth[100];
double P0_smooth_var[100];

media_blocchi (C0, C0_blocchi, 10000, 100);
media_blocchi (P0, P0_blocchi, 10000, 100);

for (int i=0; i<100; i++){
	C0_smooth[i]=media(C0_blocchi,(i+1));
	C0_smooth_var[i]=varianza(C0_blocchi,(i+1))/(i+1);
	P0_smooth[i]=media(P0_blocchi,(i+1));
	P0_smooth_var[i]=varianza(P0_blocchi,(i+1))/(i+1);
}

stampa_file("diretto_C_media",C0_smooth,100);
stampa_file("diretto_C_var",C0_smooth_var,100);
stampa_file("diretto_P_media",P0_smooth,100);
stampa_file("diretto_P_var",P0_smooth_var,100);




//SAMPLING DISCRETO

for(int i=0; i<10000; i++){
	double St=S0;
	for (int j=0; j<100; j++){
		St=St*exp((r-sigma*sigma/2)*(T/100.)+sigma*rnd.Gauss(0,1)*pow((T/100.),0.5));
	}
	ST=St;
	
	if ((ST-K)>0.){x=(ST-K);}
	else{x=0;}
	C0[i]=exp(-r*T)*x;

	if ((K-ST)>0.){x=(K-ST);}
	else{x=0;}
	P0[i]=exp(-r*T)*x;
}

cout<<"call: "<<media(C0,10000)<<"+-"<<pow(varianza(C0,10000),0.5)<<endl;
cout<<"put: "<<media(P0,10000)<<"+-"<<pow(varianza(P0,10000),0.5)<<endl;

media_blocchi (C0, C0_blocchi, 10000, 100);
media_blocchi (P0, P0_blocchi, 10000, 100);

for (int i=0; i<100; i++){
	C0_smooth[i]=media(C0_blocchi,(i+1));
	C0_smooth_var[i]=varianza(C0_blocchi,(i+1))/(i+1);
	P0_smooth[i]=media(P0_blocchi,(i+1));
	P0_smooth_var[i]=varianza(P0_blocchi,(i+1))/(i+1);
}

stampa_file("discreto_C_media",C0_smooth,100);
stampa_file("discreto_C_var",C0_smooth_var,100);
stampa_file("discreto_P_media",P0_smooth,100);
stampa_file("discreto_P_var",P0_smooth_var,100);


return 0;

}
