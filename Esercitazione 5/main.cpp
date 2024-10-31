
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
    if(argc<3){
        cerr<<"Errore: uso del programma "<<argv[0]<<" <Numero_di_blocchi> <Numero_di_step_per_blocco> <unif/gauss>"<<endl;
        return -1;
    }

    vec pos = {0.,0.,0.};
    hydro_210 distr;
    hydro_100 rtsid;
    
    string* type;
    type = new string[1];
    type[0]= argv[3];
    
    int L=atoi(argv[1]);                        //Numero di blocchi
    int O=atoi(argv[2]);                        //Numero di passi per blocco
    

    metropolis metro(distr,type[0],"../Librerie/Random Generator/");
    prova sys(L,O,1,3,type,metro,pos);
    
    sys.blocks(true);
    cout<<endl;
    
    type[0]= argv[3];
    metropolis ortem(rtsid,type[0],"../Librerie/Random Generator/");
    prova ssy(L,O,1,3,type,ortem,pos);
    
    ssy.blocks(true);
    cout<<endl;
    
    return 0;
}

