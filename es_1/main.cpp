#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "random.h"


using namespace std;

//FUNZIONI

double media (double v[], int N) {
	double sum=0;
	for (int i=0; i<N; i++){
		sum=sum+v[i];
	}
	return double(sum/N);
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

void media_blocchi (double v[], double risultati[], int N, int n) {
	if (N%n!=0){cout<<"errore Nblocchi non divisibili"<<endl;}
	int k=0;
	for (int i=0; i<(N/n); i++){
		double sum=0.;
		for (int j=0; j<n; j++){
			sum=double(sum+v[k]);
			k++;
		}
		risultati[i]=double(sum/n);
	}

}

void varianza_blocchi (double v[], double risultati[], int N, int n) {
	if (N%n!=0){cout<<"errore Nblocchi non divisibili"<<endl;}
	double mu_blocchi[N/n];
	media_blocchi(v,mu_blocchi,N,n);
	double v2[N];
	double mu_blocchi2[N/n];
	for (int i=0; i<N; i++){
		v2[i]=double(v[i]*v[i]);
	}
	media_blocchi(v2,mu_blocchi2,N,n);
	for (int i=0; i<(N/n); i++){
		risultati[i]=double(mu_blocchi2[i]-(mu_blocchi[i]*mu_blocchi[i]));
	}
}

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


//MAIN

int main (int argc, char *argv[]){


//NUMERI RANDOM
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
PARTE 1

trova:
- media e varianza a blocchi smooth
- chi quadro per più prove
per numeri pseudocasuali
*/

//MEDIA E VARIANZA

double random1[100000];
double medie_blocchi[1000];
double var_blocchi[1000];
double medie_smooth[999];
double var_smooth[999];

for (int i=0; i<100000; i++){
	random1[i]=rnd.Rannyu();
}

media_blocchi(random1,medie_blocchi,100000,100);
varianza_blocchi(random1,var_blocchi,100000,100);

for (int i=0; i<999; i++){
	medie_smooth[i]=media(medie_blocchi,i+1);
	var_smooth[i]=media(var_blocchi,i+1);
}

const char *nome_medie="medie_smooth";
const char *nome_var="var_smooth";
stampa_file(nome_medie, medie_smooth, 999);
stampa_file(nome_var, var_smooth, 999);


//CHI QUADRO

int M=10000;
double random[M];
double chi[100]={0};
ofstream chi_2;
chi_2.open("chi_2");

for (int k=0; k<100; k++){
int cont[100]={0};
	for (int i=0.0; i<M; i++){
		random[i]=rnd.Rannyu();
	}

	for (int i=0.0; i<100; i++){
		for (int j=0.0; j<M; j++){
			if (random[j]>(i/100.) && random[j]<((i+1)/100.)){
				cont[i]++;
			}
		}
		chi[k]=chi[k]+(cont[i]-100)*(cont[i]-100)/100;
	}

if (chi_2.is_open()){
	chi_2<<chi[k]<<endl;	
	}else cerr<<"problema: non si è aperto il file chi_2"<<endl;
}
chi_2.close();


/*
PARTE 2

aumentando le componenti dei blocchi, trova a che distribuzione convergono le distribuzioni:
- gaussiana
- essponenziale
- lorentziana
*/

double gauss[10000];
double exp[10000];
double lorentz[10000];
for (int i=0; i<10000; i++){
	gauss[i]=rnd.Gauss(0.,1.);
	exp[i]=rnd.Exp(1.);
	lorentz[i]=rnd.Lorentz(1.,0.);
}

double gauss_2[5000];
double exp_2[5000];
double lorentz_2[5000];
media_blocchi(gauss,gauss_2,10000,2);
media_blocchi(exp,exp_2,10000,2);
media_blocchi(lorentz,lorentz_2,10000,2);

double gauss_10[1000];
double exp_10[1000];
double lorentz_10[1000];
media_blocchi(gauss,gauss_10,10000,10);
media_blocchi(exp,exp_10,10000,10);
media_blocchi(lorentz,lorentz_10,10000,10);

double gauss_100[100];
double exp_100[100];
double lorentz_100[100];
media_blocchi(gauss,gauss_100,10000,100);
media_blocchi(exp,exp_100,10000,100);
media_blocchi(lorentz,lorentz_100,10000,100);

const char *nome;

nome="gauss_1";
stampa_file(nome, gauss, 10000);

nome="gauss_2";
stampa_file(nome, gauss_2, 5000);

nome="gauss_10";
stampa_file(nome, gauss_10, 1000);

nome="gauss_100";
stampa_file(nome, gauss_100, 100);

nome="exp_1";
stampa_file(nome, exp, 10000);

nome="exp_2";
stampa_file(nome, exp_2, 5000);

nome="exp_10";
stampa_file(nome, exp_10, 1000);

nome="exp_100";
stampa_file(nome, exp_100, 100);

nome="lorentz_1";
stampa_file(nome, lorentz, 10000);

nome="lorentz_2";
stampa_file(nome, lorentz_2, 5000);

nome="lorentz_10";
stampa_file(nome, lorentz_10, 1000);

nome="lorentz_100";
stampa_file(nome, lorentz_100, 100);

/*
PARTE 3

stima di pigreco
PI=2*L*prove/(d*successi)
- L = lunghezza ago
- d = distanza tra le rette del piano
*/

double l=7.;        
double d=10.;          
M=100000;          
int N=100;
int n=M/N;        

double vettPos[M];    
double vettL[M];      
double Nhit[M];
	

for(int i=0; i<M; i++){
	vettPos[i]=rnd.Rannyu(0.,d);
	int cont=0;
	while(cont==0){
		double x=rnd.Rannyu(0.,l/2.);
		double y=rnd.Rannyu(0.,l/2.);
		if((y*y+x*x)<(l/2.)){
			cont=1;
			vettL[i]=l/2.*(y/sqrt(pow(x,2.)+pow(y,2.)));
			}
		}
		if((vettPos[i]+vettL[i])>d || (vettPos[i]-vettL[i])<0.){
			Nhit[i]=1;

		}
		else{
			Nhit[i]=0;
		}
	}

double PI[N];
double PIsmooth[N-1];
double PIvarsmooth[N-1];

for (int i=0; i<N; i++){
	int successi=0;
	for (int j=0; j<n; j++){
		successi=successi+Nhit[i*n+j];
	}
	PI[i]=2.*l*n/(d*successi);
}

for (int i=0; i<(N-1); i++){
	PIsmooth[i]=media(PI,i+1);
	PIvarsmooth[i]=varianza(PI,i+1);	
}

stampa_file("PI_media", PIsmooth ,99);
stampa_file("PI_var", PIvarsmooth ,99);

return 0;
}



