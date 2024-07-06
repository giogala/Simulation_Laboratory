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
    if(argc<3){
            cerr<<"Errore: uso del programma "<<argv[0]<<" <input_directory> <output_directory>"<<endl;
            return -1;
        }
    //cerr<<"-2";
    int nconf = 1;
    System SYS;
    //cerr<<"-1";
    SYS.initialize(argv[1],argv[2]);
    SYS.initialize_properties();
    SYS.block_reset(0);
    //cerr<<"0";
    for(int k=0;k<10;k++){
        SYS.initialize_velocities(k);
        cerr<<"=";
        for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
            for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
                SYS.step();
                SYS.measure();
                if(j%10 == 0){
                    //          SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
                    nconf++;
                }
            }
            SYS.averages(i+1);
            SYS.block_reset(i+1);
        }
        SYS.finalize();
        SYS.read_configuration();
    }
    cerr<<endl;
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
