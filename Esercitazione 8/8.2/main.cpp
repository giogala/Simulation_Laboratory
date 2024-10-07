
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <armadillo>

#include "../../Librerie/datablocking.h"
//#include "../../Librerie/funzione.h"
//#include "../../Librerie/metropolis.h"
#include "../../Librerie/random.h"

using namespace std;
using namespace arma;


int main (int argc, char** argv) {
    if(argc<2){
        cerr<<"Errore: uso del programma "<<argv[0]<<" <input file>"<<endl;
        return -1;
    }
    
    string input = argv[1];
    
    string *type;
    type = new string[1];
    type[0]= "unif";
    //type[1]= "energy";
    
    string *prop;
    prop = new string[1];
    prop[0] = "best";
    //prop[1] = "mu";
    
    // a parametri fissati, calcolo E
    int L = int(SetProp(input,"NBLOCKS"));               //Numero di blocchi
    int O = int(SetProp(input,"NSTEPS"));                //Numero di passi per blocco
    // movimento nello spazio dei parametri
    int T = int(SetProp(input,"NTEMP"));                 //Numero di temperature esplorate
    int step = int(SetProp(input,"NTSTEP"));             //Numero di punti esplorati a T fissato
    
    double t = SetProp(input,"T0");                 //Temperatura iniziale
    double s = SetProp(input,"MOVE");               //simga di boltzmann
    double div = SetProp(input,"DIV");              //divisore della temperatura
    
    vec par = {1.,1.};                  // parametri sigma e mu inziali
    vec pos = {0.2};                    // posizione inziale di campionamento della funzione d'onda
    hat distr(1.,1.);                   // funzione d'onda inizializzata con sigma = mu = 1.
    
    metropolis metro(distr,type[0],"../../Librerie/Random Generator/");
    energy sys(L,O,1,1,type,metro,pos);
    boltzmann boltz(t,sys);
    metropolis bigM(boltz,type[0],"../../Librerie/Random Generator/");
    bigM.SetS(s);
    annealing null(T,step,1,2,div,prop,sys,bigM);
    
    null.blocks(true);
    cout<<endl;

    return 0;
}

