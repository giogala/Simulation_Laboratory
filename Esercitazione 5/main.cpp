
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <armadillo>

#include "funzione.h"
#include "metropolis.h"
#include "random.h"

using namespace std;
using namespace arma;


int main (int argc, char** argv) {
    if(argc<3){
        cerr<<"Errore: uso del programma "<<argv[0]<<" <Numero_di_passi> <Numero_blocchi> <distribuzione_di_transizione>"<<endl;
        return -1;
    }
    vec p_i = {0.,0.,0.};
    vec pos = {0.,0.,0.};
    hydro_210 distr;
    string type = argv[3];
    //if(type == "gauss") gauss3d trans;
    //if(type == "unif") unif trans;
    //funzione& trans;
    metropolis metro(distr,type,p_i);
    //File di output
    ofstream fout;
    fout.open("Data.txt");

    int M=atoi(argv[1]);                        //Numero di passi
    int L=atoi(argv[2]);                        //Numero di blocchi
    int O=M/L;                                  //Numero di repliche per blocco
    
    double s=0; double s2=0; double t=0; double t2=0;
    // iterazione sui blocchi in cui ho diviso le repliche
    for(int k=0;k<L;k++){
        //iterazione sulle repliche all'interno di un blocco
        for(int j=0;j<O;j++){
            pos = metro.move(pos);
            s+= sqrt(dot(pos,pos))/O;
            }
        metro.adjust();
        s2=s*s;
        t=(t*k+s)/(k+1);
        t2=(t2*k+s2)/(k+1);
        s2=0; s=0;
        double e = sqrt((t2-t*t)/k);      // il "/L" è dovuto al fatto che sto mediando su L blocchi
        if(k==0) e=0;
        fout<<t<<"\t"<<e<<endl;
    }

    
    return 0;
}

