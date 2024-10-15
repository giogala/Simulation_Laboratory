
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
    int f = SetProp("input.txt","GEN");
    string t = SetType("input.txt","TYPE");
    Random ran;
    ran.initRnd("../../Librerie/Random Generator/");
    population test(N,ran,"../"+t+"/");
    if(t=="CIRCLE") test.Circle(p);
    if(t=="SQUARE") test.Square(p);
    test.Check();
    
    for(int k=0;k<f;k++){
        test.Mutate();
        test.Evolve(50);
        
        Progress_Bar(k,f);
        test.Check();
        test.Sort();
        test.Print(k);
        test.L2(k);
    }
    
    return 0;
}

