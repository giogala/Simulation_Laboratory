#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "../../Librerie/random.h"
#include "../../Librerie/library.h"

using namespace std;

int main (int argc, char** argv) {
    if(argc<2){
        cerr<<"Errore: uso del programma "<<argv[0]<<" <Numero_di_lanci>"<<endl;
        return -1;
    }
    
    Random rnd;                                 //Imposto il generatore di numeri casuali
    rnd.initRnd("../../Librerie/Random Generator/");
    
    int N=atoi(argv[1]);
    int M[4]={1,2,10,100};
    
    //File di output - uno per tipo di distribuzione
    ofstream fout("Data.txt");
    ofstream u_out("Data_u.txt");
    ofstream e_out("Data_e.txt");
    ofstream c_out("Data_c.txt");
    
    //Parte 1
    for(int i=0;i<N;i++){
        double a=rnd.Exp(1.);
        double b=rnd.CauLor(0.,1.);
        fout<<a<<"\t"<<b<<"\t"<<a+b<<endl;
        Progress_Bar(i,N);
    }
    //Parte 2
    for(int i=0;i<N;i++){
        u_out<<i;
        e_out<<i;
        c_out<<i;
        for(int j=0;j<4;j++){
            double unif=0;
            double exp=0;
            double cau=0;
            for(int k=0;k<M[j];k++){                    //genero numeri con ciascuna distribuzione
                unif+=rnd.Rannyu()/M[j];                //ne sommo 1,2,10 o 100 a seconda del ciclo
                exp+=rnd.Exp(1.)/M[j];
                cau+=rnd.CauLor(0.,1.)/M[j];
            }
            u_out<<"\t"<<unif;                          //stampo su file (divisi per distribuzione)
            e_out<<"\t"<<exp;                           //ciascuna colonna è data da somma di 2, 10 ...
            c_out<<"\t"<<cau;
        }
        u_out<<endl;
        e_out<<endl;
        c_out<<endl;
        Progress_Bar(i,N);
    }
    
    fout.close();
    u_out.close();
    e_out.close();
    c_out.close();
    return 0;
}

