
#define _USE_MATH_DEFINES
#include "mpi.h"
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
M_PI

int main (int argc, char** argv) {
    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int p = SetProp("input.txt","NPOINTS");
    int f = SetProp("input.txt","GEN");
    int a = SetProp("input.txt","MOVE");
    int m = SetProp("input.txt","MIGR");
    string t = SetType("input.txt","TYPE");
    
    Random ran;
    ran.initRnd("../../Librerie/Random Generator/");
    population test("input.txt",ran,"../"+t+"/");
    
    if(t=="CIRCLE") test.Circle(p);
    else if(t=="SQUARE") test.Square(p);
    else if(t=="FILE") {
        test.Rnd().initParallel("../../Librerie/Random Generator/",rank);
        test.File(SetType("input.txt","FILE"));
    }
    else{
        cerr<<"Unknown initial configuration"<<endl;
        return -1;
    }
    test.Rnd().initParallel("../../Librerie/Random Generator/",rank);
    test.Check();
    
    for(int k=0;k<f;k++){
        test.Mutate();
        test.Evolve();
        test.Select();
        
        if(rank == 0)Progress_Bar(k,f);
        test.Check();
        test.Sort();
        if(rank == 0)test.Print(k);
        if(rank == 0)test.L2(k);
        if(k%a == 0) test.Migration(rank,size,m);
    }
    
    MPI_Finalize();
    return 0;
}

