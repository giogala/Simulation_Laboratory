
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <armadillo>

#include "../../Librerie/population.h"
#include "../../Librerie/datablocking.h"
#include "../../Librerie/funzione.h"
#include "../../Librerie/metropolis.h"
#include "../../Librerie/random.h"

using namespace std;
using namespace arma;


int main (int argc, char** argv) {
    if(argc<2){
        cerr<<"Errore: uso del programma "<<argv[0]<<" <circle/square>"<<endl;
        return -1;
    }
    int N = SetProp("input.txt","NELEM");                   //Numero di elementi
    int p = SetProp("input.txt","NPOINTS");
    int i = SetProp("input.txt","PRINT");
    Random ran;
    ran.initRnd("../../Librerie/Random Generator/");
    population test(N,ran);
    //if(argv[1]=="circle") test.Circle(p);
    test.Circle(p);
    test.Xover(0,i);
    test.Check();
    
    test.Sort();
    test.L2();
    test.Print();
    
    return 0;
}

