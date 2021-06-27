#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <ostream>
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

double phi (double x, double mu, double sigma){
	//double phi=exp(-pow(x-mu,2)/(sigma*sigma))+exp(-(x*x+mu*mu)/(sigma*sigma))+exp(-pow(x+mu,2)/(sigma*sigma));
	double phi=exp(-pow(x-mu,2)/(2*sigma*sigma))+exp(-pow(x+mu,2)/(2*sigma*sigma));
	return phi;
}

double Tphi (double x, double mu, double sigma){
	double Tphi=(exp(-pow(x-mu,2)/(2*sigma*sigma))*(1-pow(x-mu,2)/(2*sigma*sigma))+exp(-pow(x+mu,2)/(2*sigma*sigma))*(1-pow(x+mu,2)/(2*sigma*sigma)))/(2*sigma*sigma);
	return Tphi;
}

double V (double x){
	return (pow(x,4)-5*pow(x,2)/2);
}

int main (){

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
   
//8.1 calcolo <H>

double l=4.4;
double mu=1.5;
double sigma=1.;
double x=0.;
double xnew;
double p;
int n_blocchi=1000;
int n_componenti=200;
double fx;
double sum=0.;
double fx_ave[n_blocchi];
double fx_ave_smooth[n_blocchi];
double fx_var_smooth[n_blocchi];
ofstream prx;

int k=0;
prx.open("psi2_x");

/*for (int i=0; i<step/10; i++){
	xnew=x+(rnd.Rannyu()*2*l-l);
	p=pow(phi(xnew,mu,sigma),2)/pow(phi(x,mu,sigma),2);
	if (p>1){p=1;}
	double rand=rnd.Rannyu();
	if(rand<p){
		x=xnew;
	}
}*/

for (int i=0; i<n_blocchi; i++){
	sum=0;
	for(int j=0; j<n_componenti; j++){
		xnew=x+(rnd.Rannyu()*2*l-l);
		p=pow(phi(xnew,mu,sigma),2)/pow(phi(x,mu,sigma),2);
		if (p>1){p=1;}
		double rand=rnd.Rannyu();
		if(rand<p){
			x=xnew;
			k++;
		}
		fx=Tphi(x,mu,sigma)/phi(x,mu,sigma)+V(x);
		sum=sum+fx;
		prx<<x<<endl;
		
	}
	fx_ave[i]=sum/n_componenti;
}
prx.close();

cout<<"<H>   ->   mu = "<<mu<<"   sigma = "<<sigma<<endl;
cout<<"probabilità accettazzione: "<<k*100/n_blocchi/n_componenti<<"%"<<endl;
cout<<"media: "<<media(fx_ave, n_blocchi)<<endl;
cout<<"err: "<<sqrt(varianza(fx_ave, n_blocchi)/n_blocchi)<<endl<<endl;

for (int i=0; i<n_blocchi; i++){
	fx_ave_smooth[i]=media(fx_ave, i+1);
	fx_var_smooth[i]=sqrt(varianza(fx_ave, i+1)/(i+1));
}

stampa_file("H_smooth",fx_ave_smooth, n_blocchi);
stampa_file("H_err_smooth",fx_var_smooth, n_blocchi);

//8.2 ricerca di mu e sigma
ofstream pr;
pr.open("H_tentativi");

double H_T, H_err_T;
for (mu=0; mu<4.1; mu+=0.1){
	for (sigma=0.1; sigma<3.1; sigma+=0.1){
		k=0;
		for (int i=0; i<n_blocchi; i++){
			sum=0;
			for(int j=0; j<n_componenti; j++){
				xnew=x+(rnd.Rannyu()*2*l-l);
				p=pow(phi(xnew,mu,sigma),2)/pow(phi(x,mu,sigma),2);
				if (p>1){p=1;}
				double rand=rnd.Rannyu();
				if(rand<p){
					x=xnew;
					k++;
				}
				fx=Tphi(x,mu,sigma)/phi(x,mu,sigma)+V(x);
				sum=sum+fx;
			}
			fx_ave[i]=sum/n_componenti;
		}
		H_T=media(fx_ave, n_blocchi);
		H_err_T=sqrt(varianza(fx_ave, n_blocchi)/n_blocchi);
		cout<<"<H>   ->   mu = "<<mu<<"   sigma = "<<sigma<<endl;
		cout<<"probabilità accettazzione: "<<k*100/n_blocchi/n_componenti<<"%"<<endl;
		cout<<"media: "<<H_T<<"   err: "<<H_err_T<<endl<<endl;
		pr<<setw(16)<<H_T<<setw(16)<<H_err_T<<setw(16)<<mu<<setw(16)<<sigma<<endl;
	}
}
pr.close();



ofstream pr2;
pr2.open("H_tentativi2");
//double H_T, H_err_T;
l=2.;
for (mu=0.; mu<0.25; mu+=0.005){
	for (sigma=0.95; sigma<1.1; sigma+=0.005){
		k=0;
		for (int i=0; i<n_blocchi; i++){
			sum=0;
			for(int j=0; j<n_componenti; j++){
				xnew=x+(rnd.Rannyu()*2*l-l);
				p=pow(phi(xnew,mu,sigma),2)/pow(phi(x,mu,sigma),2);
				if (p>1){p=1;}
				double rand=rnd.Rannyu();
				if(rand<p){
					x=xnew;
					k++;
				}
				fx=Tphi(x,mu,sigma)/phi(x,mu,sigma)+V(x);
				sum=sum+fx;
			}
			fx_ave[i]=sum/n_componenti;
		}
		H_T=media(fx_ave, n_blocchi);
		H_err_T=sqrt(varianza(fx_ave, n_blocchi)/n_blocchi);
		cout<<"<H>   ->   mu = "<<mu<<"   sigma = "<<sigma<<endl;
		cout<<"probabilità accettazzione: "<<k*100/n_blocchi/n_componenti<<"%"<<endl;
		cout<<"media: "<<H_T<<"   err: "<<H_err_T<<endl<<endl;
		pr2<<setw(16)<<H_T<<setw(16)<<H_err_T<<setw(16)<<mu<<setw(16)<<sigma<<endl;
	}
}
pr2.close();


//PROVA CON TERMALIZZAZIONE
/*mu=rnd.Rannyu(0.,2.);
sigma=rnd.Rannyu(0.,2.5);
l=2.;
x=0;
for (int i=0; i<n_blocchi; i++){
	sum=0;
	for(int j=0; j<n_componenti; j++){
		xnew=x+(rnd.Rannyu()*2*l-l);
		p=pow(phi(xnew,mu,sigma),2)/pow(phi(x,mu,sigma),2);
		if (p>1){p=1;}
		double rand=rnd.Rannyu();
		if(rand<p){x=xnew;}
		fx=Tphi(x,mu,sigma)/phi(x,mu,sigma)+V(x);
		sum=sum+fx;	
	}
	fx_ave[i]=sum/n_componenti;
}
double H= media(fx_ave, n_blocchi);
double newmu, newsigma, newH;
double bestmu=mu;
double bestsigma=sigma;
double bestH=H;

for (int T=0; T<1; T++){
	cout<<T<<" di 100"<<endl;
	for (int k=0; k<500; k++){
		newmu=mu+rnd.Rannyu((-0.5)/(T+1), (0.5)/(T+1));
		newsigma=sigma+rnd.Rannyu((-0.5)/(T+1), (0.5)/(T+1));
		for (int i=0; i<n_blocchi; i++){
			sum=0;
			for(int j=0; j<n_componenti; j++){
				xnew=x+(rnd.Rannyu()*2*l-l);
				p=pow(phi(xnew,mu,sigma),2)/pow(phi(x,mu,sigma),2);
				if (p>1){p=1;}
				double rand=rnd.Rannyu();
				if(rand<p){x=xnew;}
				fx=Tphi(x,mu,sigma)/phi(x,mu,sigma)+V(x);
				sum=sum+fx;	
			}	
			fx_ave[i]=sum/n_componenti;
		}
		newH= media(fx_ave, n_blocchi);
		double p2=exp(-T*(newH-H));
		if(newH<H){p=1;}
		if(rnd.Rannyu()<p2){
			mu=newmu;
			sigma=newsigma;
			H=newH;
		}
		if(newH<bestH){
			bestmu=newmu;
			bestsigma=newsigma;
			bestH=newH;
			cout<<"best H: "<<bestH<<"  "<<mu<<" "<<sigma<<endl;
		}
	}
}
*/

//valor medio <H>
mu=0.09;
sigma=1.03;
l=2.;
x=0.;
sum=0.;
k=0;
ofstream prx2;
prx2.open("psi2_x_GS");
for (int i=0; i<n_blocchi; i++){
	sum=0;
	for(int j=0; j<n_componenti; j++){
		xnew=x+(rnd.Rannyu()*2*l-l);
		p=pow(phi(xnew,mu,sigma),2)/pow(phi(x,mu,sigma),2);
		if (p>1){p=1;}
		double rand=rnd.Rannyu();
		if(rand<p){
			x=xnew;
			k++;
		}
		fx=Tphi(x,mu,sigma)/phi(x,mu,sigma)+V(x);
		sum=sum+fx;
		prx2<<x<<endl;
		
	}
	fx_ave[i]=sum/n_componenti;
}
prx2.close();

cout<<"<H>   ->   mu = "<<mu<<"   sigma = "<<sigma<<endl;
cout<<"probabilità accettazzione: "<<k*100/n_blocchi/n_componenti<<"%"<<endl;
cout<<"media: "<<media(fx_ave, n_blocchi)<<endl;
cout<<"err: "<<sqrt(varianza(fx_ave, n_blocchi)/n_blocchi)<<endl<<endl;

for (int i=0; i<n_blocchi; i++){
	fx_ave_smooth[i]=media(fx_ave, i+1);
	fx_var_smooth[i]=sqrt(varianza(fx_ave, i+1)/(i+1));
}

stampa_file("H_smooth_GS",fx_ave_smooth, n_blocchi);
stampa_file("H_err_smooth_GS",fx_var_smooth, n_blocchi);

rnd.SaveSeed();
return 0;

}

