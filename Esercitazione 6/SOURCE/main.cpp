/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
// equilibration main for T in 0.5 1.5
#include <iostream>
#include "system.h"

using namespace std;

int main (int argc, char *argv[]){
    if(argc<1){
        cerr<<"Error! Usage of the program: "<<argv[0]<<" <input.dat>"<<endl;
        return -1;
    }
    int nconf = 1;
    System SYS;
    SYS.initialize("../INPUT","../OUTPUT","../../Librerie/Random Generator",argv[1]);
    SYS.initialize_properties();
    SYS.block_reset(0);
    SYS.read_configuration();
    //SYS.initialize_velocities(0);

    for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
        for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
            SYS.step();
            SYS.measure();
            if(j%10 == 0){
                //SYS.write_XYZ(nconf);
                //Write actual configuration in XYZformat //Commented to avoid "filesystem full"!
                nconf++;
            }
        }
        SYS.averages(i+1);
        SYS.block_reset(i+1);
    }
    SYS.finalize();
    //SYS.read_configuration();
    
    return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
