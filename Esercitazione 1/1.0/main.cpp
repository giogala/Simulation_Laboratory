#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "../../Librerie/Random Generator/random.h"

using namespace std;

int main (int argc, char** argv) {
    if(argc<3){
        cerr<<"Errore: uso del programma "<<argv[0]<<" <Numero_di_lanci> <Numero_di_blocchi>"<<endl;
        return -1;
    }
    
    Random rnd;                                 //Imposto il generatore di numeri casuali
    int seed[4]={1,2,3,4};
    //rnd.SetRandom("../../Librerie/Random Generator/Primes","../../Librerie/Random Generator/seed.in");
    rnd.SetRandom(seed,5,6);

    
    int N=atoi(argv[1]);                        //Numero di lanci casuali
    int M=atoi(argv[2]);                        //Numero di blocchi in cui dividere il dataset
    int L=N/M;                                  //Numero di lanci per blocco
    
    double c=0;                                 //contatore per i valori medi del blocco
    double c2=0;                                //idem ma ^2
    double t=0;                                 //contatore la somma dei valori medi fino all'i-esimo blocco
    double t2=0;                                //idem ma ^2

    ofstream fout("Data.txt");                  //File di output
    
    for(int i=0;i<M+1;i++){
        for(int j=0;j<L+1;j++){
            c+= rnd.Rannyu()/L;                 //Genero il numero casuale
        }
        c2=pow(c,2);
        t=(t*i+c)/(i+1);                        //preso il valor medio dei primi i blocchi, lo motiplico per i,
        t2=(t2*i+c2)/(i+1);                     //aggiungo l'(i+1)-esimo valore e quindi divido per (i+1)
                                                //stessa cosa per il valor medio dei <r^2>
        double err=0;
        if(i!=0){
            err= sqrt((-pow(t,2)+t2)/i);        //Calcolo dell'errore
        }
        fout<<t<<"\t"<<err<<endl;
        c=0; c2=0;                              //resetto i contatori per passare al blocco successivo
    }
    fout.close();
    
    return 0;
}

