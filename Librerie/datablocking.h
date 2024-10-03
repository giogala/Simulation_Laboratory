//
//  datablocking.h
//  
//
//  Created by Giovanni Galafassi on 02/10/24.
//
#define _USE_MATH_DEFINES
#ifndef __datablocking_h_
#define __datablocking_h_
#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdio.h>
#include "metropolis.h"
#include "library.h"
#include "funzione.h"
#include "random.h"

using namespace std;
using namespace arma;

class datablocking{
public:
    datablocking(int L, int O, int p, int dim,string* file){
        _nblock = L;
        _nstep = O;
        _nprop = p;
        _f = file;
        
        _dim = dim;
        _state.resize(dim);
        _state.zeros();
        
        _av.resize(p);
        blk_av.resize(p);
        glb_av.resize(p);
        glb_2.resize(p);
        _av.zeros();
        blk_av.zeros();
        glb_av.zeros();
        glb_2.zeros();
    };
    virtual ~datablocking(){;};
    
    void blocks(){
        for(int i=0;i<_nblock;i++){
            blk_av.zeros();
            _state = increase(_state);
            Progress_Bar(i+1,_nblock);
            for(int j=0;j<_nprop;j++) {
                averages(j);
                print(i+1,j);
            }
        }
    };
    
    void averages(int prop){
        _av(prop) = blk_av(prop)/double(_nstep);
        adjust(prop);
        glb_av(prop) += _av(prop);
        glb_2(prop) += _av(prop) * _av(prop);
    };
    
    void print(int blk,int prop){
        ofstream fout;
        fout.open(_f[prop]+".txt",ios::app);
        average  = _av(prop);
        sum_average = glb_av(prop);
        sum_ave2 = glb_2(prop);
        fout << blk
          << "\t" << average
          << "\t" << sum_average/double(blk)
          << "\t" << this->error(sum_average, sum_ave2, blk) << endl;
        fout.close();
    };
    
    double error(double acc, double acc2, int blk){
        if(blk <= 1) return 0.0;
        else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
    };

    virtual vec increase(vec pos) =0;
    virtual void adjust(int prop) =0;
    
protected:
    int _nblock;
    int _nstep;
    vec _state;
    int _dim;
    vec _av;
    vec blk_av;
    vec glb_av;
    vec glb_2;
    double average, sum_average, sum_ave2;
    int _nprop;
    string* _f;
};

class prova: public datablocking{
public:
    prova(int L, int O, int p, int dim,string* file, metropolis& metro,vec pos):datablocking(L,O,p,dim,file),_metro(metro){
        string add;
        if(_metro.GetDistr().whoamI()=="hydro_100") add = "100";
        else if(_metro.GetDistr().whoamI()=="hydro_210") add = "210";
        else{
            cerr<<"Attenzione! Distribuzione non riconosciuta, pericolo sovrascrizione"<<endl;
        }
        for(int i=0;i<_nprop;i++){
            ofstream fout;
            _f[i]=_f[i]+"_"+add;
            fout.open(_f[i]+".txt");
            fout<<"BLOCK\tACTUAL_P\tP_AVE\tERROR"<<endl;
            fout.close();
        }
        _state = pos;
    };
    vec increase(vec pos){
        for(int k=0;k<_nstep;k++){
            pos = _metro.move(pos);
            blk_av(0) += sqrt(dot(pos,pos));
        }
        return pos;
    };
    void adjust(int prop){
        if(prop==0)_metro.adjust();
    };
private:
    metropolis& _metro;
};

class energy: public datablocking{
public:
    energy(int L, int O, int p, int dim,string* file, metropolis& metro,vec pos):datablocking(L,O,p,dim,file),_metro(metro){
        for(int i=0;i<_nprop;i++){
            ofstream fout;
            fout.open(_f[i]+".txt");
            fout<<"BLOCK\tACTUAL_P\tP_AVE\tERROR"<<endl;
            fout.close();
        }
        _state = pos;
    };
    vec increase(vec pos){
        for(int k=0;k<_nstep;k++){
            pos = _metro.move(pos);
            blk_av(0) += _metro.GetDistr().der(pos);
        }
        return pos;
    };
    void adjust(int prop){
        if(prop==0)_metro.adjust();
    };
private:
    metropolis& _metro;
};



#endif /* datablocking_h */
