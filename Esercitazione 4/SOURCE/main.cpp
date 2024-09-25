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
#include "../../Librerie/system.h"
#include "../../Librerie/library.h"

using namespace std;

int main (int argc, char *argv[]){
    if(argc<3){
        cerr<<"Error! Usage of the program: "<<argv[0]<<" <input.dat> <y/n> (if you want or not XYZ format) "<<endl;
        return -1;
    }
    int nconf = 1;
    bool xyz = false;
    if(strcmp(argv[2],"y") == 0) xyz = true;
    System SYS;
    SYS.initialize("../INPUT","../OUTPUT","../../Librerie/Random Generator",argv[1]);
    SYS.initialize_properties();
    SYS.block_reset(0);
    SYS.initialize_velocities();
    for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
        for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
            SYS.step();
            SYS.measure();
            if(j%100 == 0 and xyz){
                SYS.write_XYZ(nconf);   //Write actual configuration in XYZformat
                                        //Commented to avoid "filesystem full"!
                nconf++;
            }
        }
        Progress_Bar(i+1,SYS.get_nbl());
        SYS.averages(i+1);
        SYS.block_reset(i+1);
    }
    SYS.finalize();
    //SYS.read_configuration();
    cout<<endl;
    
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
