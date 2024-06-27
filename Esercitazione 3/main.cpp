
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "../Librerie/Random Generator/random.h"
#include "../Librerie/Random Generator/posizione.h"

using namespace std;

double max(double a, double b){
    if(a>b) return a;
    else return b;
}

int main (int argc, char** argv) {
    if(argc<3){
        cerr<<"Errore: uso del programma "<<argv[0]<<" <Numero_di_passi> <Numero_di_repliche> <Numero_blocchi_di_repliche>"<<endl;
        return -1;
    }
    
    Random rnd;                                 //Imposto il generatore di numeri casuali
    int seed[4]={5,2,3,4};
    rnd.SetRandom(seed,5,100);
    
    //File di output
    ofstream fout[4];
    for(int i=0;i<4;i++){
        fout[i].open("Data"+to_string(i)+".txt");
        fout[i]<<"Price"<<"\t"<<"Error"<<endl;
    }
    
    int N=atoi(argv[1]);                        //Numero di passi
    int M=atoi(argv[2]);                        //Numero di repliche
    int L=atoi(argv[3]);                        //Numero di blocchi in cui suddividerwe le M repliche
    int O=M/L;                                  //Numero di repliche per blocco
    
    // Vettore con tutte le posizioni di tutte le repliche del RW inizializzate nell'origine
    //vector <vector <double> > pos;
    double s0=100.;
    double K=100.;
    double r=0.1;
    double v=0.25;
    double T=1.;
    double s[4]={0,0,0,0};             //variabile locale nel blocco per il calcolo del valor medio
    double s2[4]={0,0,0,0};            //variabile locale nel blocco per il calcolo della varianza
    double t[4]={0,0,0,0};             //variabile globale per il calcolo del valor medio dei valori medi a i fissato
    double t2[4]={0,0,0,0};            //variabile globale per il calcolo della dev. std della media a i fissato
    
    
    //iterazione sulla modalità di calcolo
    for(int i=0;i<4;i++){
        // iterazione sui blocchi in cui ho diviso le repliche
        for(int k=0;k<L;k++){
            //iterazione sulle repliche all'interno di un blocco
            for(int j=0;j<O;j++){
                if(i==0) s[0]+=exp(-r*T)*max(0,s0*exp((r-v*v/2)*T+v*rnd.Gauss(0.,1.)*sqrt(T))-K)/O;
                if(i==1) s[1]+=exp(-r*T)*max(0,K-s0*exp((r-v*v/2)*T+v*rnd.Gauss(0.,1.)*sqrt(T)))/O;
                else{
                    double x=s0;
                    for(int l=0;l<N;l++){
                        x+=(double)(r*T/N+v*rnd.Gauss(0.,1.)*sqrt(T/N));
                    }
                    if(i==2)s[2]+=exp(-r*T)*max(0,x-K);
                    if(i==3)s[3]+=exp(-r*T)*max(0,K-x);
                }
            }
            s2[i]=s[i]*s[i];
            t[i]=(t[i]*k+s[i])/(k+1);
            t2[i]=(t2[i]*k+s2[i])/(k+1);
            s2[i]=0; s[i]=0;
            
            double e = sqrt((t2[i]-t[i]*t[i])/k);      // il "/L" è dovuto al fatto che sto mediando su L blocchi
            if(k==0) e=0;
            fout[i]<<t[i]<<"\t"<<e<<endl;
        }
        fout[i].close();
    }
    return 0;
}

