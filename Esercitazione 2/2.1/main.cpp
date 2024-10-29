
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "../../Librerie/random.h"
#include "../../Librerie/library.h"

using namespace std;

int main (int argc, char** argv) {
    if(argc<2){
        cerr<<"Errore: uso del programma "<<argv[0]<<" <Numero_di_lanci> <Numero_di_Blocchi"<<endl;
        return -1;
    }
    
    Random rnd;                                 //Imposto il generatore di numeri casuali
    rnd.initRnd("../../Librerie/Random Generator/");
    
    //File di output
    ofstream fout("Data.txt");
    ofstream fout2("Data2.txt");

    
    int N=atoi(argv[1]);                        //Numero di lanci casuali
    int M=atoi(argv[2]);                        //Numero di blocchi in cui dividere il dataset e cicli di p3
    int L=N/M;                                  //Numero di lanci per blocco
    
    //contatore integrale sul singolo blocco
    double I=0;     //uniforme
    double I2=0;
    double q=0;     //importance sampling
    double q2=0;
    //contatore integrale globale
    double mean=0;  //uniforme
    double var=0;
    double nem=0;   //importance sampling
    double rav=0;
    
    for(int j=0;j<M+1;j++){
        for(int i=0;i<L;i++){
            I+=(M_PI/2)*cos(M_PI*rnd.Rannyu()/2)/L;
            double y = 1-sqrt(1-rnd.Rannyu());      //generatore di numeri secondo la p(x) = 2(1-x)
            q+= M_PI*cos(M_PI*y/2)/((1-y)*4*L);     //sommo L valori casuali della g(x) secondo la p(x)
        }
        I2=pow(I,2);
        q2=pow(q,2);
        
        mean=(mean*j+I)/(j+1);
        var=(var*j+I2)/(j+1);
        nem=(nem*j+q)/(j+1);
        rav=(rav*j+q2)/(j+1);
        
        double e1 = sqrt((var-mean*mean)/j);
        double e2 = sqrt((rav-nem*nem)/j);
        if(j==0){
            e1 = 0.; e2 = 0.;
        }
        
        fout<<j<<"\t"<<mean<<"\t"<<e1<<endl;
        fout2<<j<<"\t"<<nem<<"\t"<<e2<<endl;
        I=0; I2=0; q=0; q2=0;
        Progress_Bar(j,M);
    }
    
    fout.close();
    fout2.close();
    return 0;
}

