#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "../../Librerie/random.h"
#include "../../Librerie/library.h"

using namespace std;

int main (int argc, char** argv) {
    if(argc<4){
        cerr<<"Errore: uso del programma "<<argv[0]<<" <Numero_di_lanci> <Numero_di_suddivisioni> <Numero_di_ripetizioni_parte3>"<<endl;
        return -1;
    }
    
    Random rnd;                                 //Imposto il generatore di numeri casuali
    rnd.initRnd("../../Librerie/Random Generator/");

    
    int N=atoi(argv[1]);                        //Numero di lanci casuali
    int M=atoi(argv[2]);                        //Numero di blocchi in cui dividere il dataset e cicli di p3
    int T=atoi(argv[3]);                        //Numero di lanci per iterazione della parte 3
    int L=N/M;                                  //Numero di lanci per blocco
    
    //contatore per i valori medi del blocco e i rispettivi quadtrati
    double d=0;             //parte 1
    double d2=0;
    double c=0;             //parte 2
    double c2=0;
    vector <int> n;         //parte 3 - contatore numeri ricaduti nel k-esimo intervallo
    
    //contatore la somma dei valori medi fino all'i-esimo blocco e rispettivi quadtrati
    double t=0;             //parte 1
    double t2=0;
    double e=0;             //parte 2
    double e2=0;
    
    //File di output
    ofstream fout("Data.txt");
    ofstream fout2("Data2.txt");
    ofstream fout3("Data3.txt");
    
    for(int i=0;i<M+1;i++){
        //Genero il numero casuale
        for(int j=0;j<L+1;j++){
            d+= rnd.Rannyu()/L;                 //parte 1
            c+= pow(rnd.Rannyu()-0.5,2)/L;      //parte 2
        }
        c2=pow(c,2);
        d2=pow(d,2);
        
        //preso il valor medio dei primi "i" blocchi, lo motiplico per i,
        //aggiungo l'(i+1)-esimo valore e quindi divido per (i+1) e così iterativamente
        //stessa cosa per il valor medio dei <r^2>
        t=(t*i+d)/(i+1);
        t2=(t2*i+d2)/(i+1);
        e=(e*i+c)/(i+1);
        e2=(e2*i+c2)/(i+1);
        
        double err=0;
        double ers=0;
        if(i!=0){
            err= sqrt((-pow(t,2)+t2)/i);        //Calcolo dell'errore
            ers= sqrt((-pow(e,2)+e2)/i);
        }
        fout<<i<<"\t"<<t<<"\t"<<err<<endl;
        fout2<<i<<"\t"<<e<<"\t"<<ers<<endl;
        c=0; c2=0; d=0; d2=0;                   //resetto i contatori per passare al blocco successivo
        
        //parte 3
        double b=0;
        double chi2i=0;
        for(int j=0;j<T;j++){
            b=rnd.Rannyu();
            for(int k=0;k<100;k++){
                if(i==0) n.push_back(0);
                if((double)k/100<=b and b<(double)(k+1)/100) n[k]++;
            }
        }
        for(int k=0;k<100;k++){
            chi2i+=pow(n[k]-T/100,2)*100/T;
            n[k]=0;
        }
        fout3<<i<<"\t"<<chi2i<<endl;
        Progress_Bar(i,M);
    }
    fout.close();
    fout2.close();
    fout3.close();
    return 0;
}

