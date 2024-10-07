
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <armadillo>

#include "../../Librerie/datablocking.h"
#include "../../Librerie/funzione.h"
#include "../../Librerie/metropolis.h"
#include "../../Librerie/random.h"

using namespace std;
using namespace arma;


int main (int argc, char** argv) {
    if(argc<2){
        cerr<<"Errore: uso del programma "<<argv[0]<<" <sigma> <mu>"<<endl;
        return -1;
    }
    vec pos = {0.2};
    hat distr(atof(argv[1]),atof(argv[2]));
    string *type;
    type = new string[1];
    type[0] = "unif";
    int L = SetProp("input.txt","NBLOCKS");                 //Numero di blocchi
    int O = SetProp("input.txt","NSTEPS");                  //Numero di passi per blocco
    
    metropolis metro(distr,type[0],"../../Librerie/Random Generator/");
    energy sys(L,O,1,1,type,metro,pos);
    
    sys.blocks(true);
    cout<<endl;


    
    return 0;
}

