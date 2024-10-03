
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <armadillo>

#include "../Librerie/datablocking.h"
#include "../Librerie/funzione.h"
#include "../Librerie/metropolis.h"
#include "../Librerie/random.h"

using namespace std;
using namespace arma;


int main (int argc, char** argv) {
    if(argc<2){
        cerr<<"Errore: uso del programma "<<argv[0]<<" <Numero_di_blocchi> <Numero_di_step_per_blocco>"<<endl;
        return -1;
    }
    vec pos = {0.};
    hat distr(1.,1.);
    string *type;
    type = new string[1];
    type[0]= "unif";
    int L=atoi(argv[1]);                        //Numero di blocchi
    int O=atoi(argv[2]);                        //Numero di passi per blocco
    
    metropolis metro(distr,type[0],"../Librerie/Random Generator/");
    energy sys(L,O,1,1,type,metro,pos);
    
    sys.blocks();
    cout<<endl;


    
    return 0;
}

