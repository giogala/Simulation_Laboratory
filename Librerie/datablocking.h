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
    
    void blocks(bool prt);
    void averages(int prop);
    virtual double print(int blk,int prop,bool prt);
    double error(double acc, double acc2, int blk);
    void reset();
    int N_blk();

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
            fout<<"BLOCK\tACTUAL_V\tV_AVE\tERROR"<<endl;
            fout.close();
        }
        _state = pos;
    };
    vec increase(vec pos);
    void adjust(int prop);
private:
    metropolis& _metro;
};

class energy: public datablocking{
public:
    energy(int L, int O, int p, int dim,string* file, metropolis& metro,vec pos):datablocking(L,O,p,dim,file),_metro(metro){
        for(int i=0;i<_nprop;i++){
            ofstream fout;
            fout.open(_f[i]+".txt");
            fout<<"BLOCK\tACTUAL_V\tV_AVE\tERROR"<<endl;
            fout.close();
        }
        _state = pos;
    };
    vec increase(vec pos);
    void adjust(int prop);
    vec pos();
    metropolis& GetMetro();
    
private:
    metropolis& _metro;
};
    
    

class annealing: public datablocking{
public:
    annealing(int L, int O, int p, int dim,double div,string* file, energy& ener,metropolis& metro):datablocking(L,O,p,dim,file),_ener(ener),_metro(metro){
        for(int i=0;i<_nprop;i++){
            ofstream fout;
            fout.open(_f[i]+".txt");
            if(i==0) fout<<"TEMP\tSIGMA\tMU\tENERGY"<<endl;
            else fout<<"BLOCK\tACTUAL_V\tV_AVE\tERROR"<<endl;
            fout.close();
        }
        _div = div;
        _t = _metro.GetDistr().GetParams(0);
        
        for (int i=0;i<_dim;i++) _state(i) = _ener.GetMetro().GetDistr().GetParams(i);
        _min.resize(dim);
        _min = _state;
        _ener.blocks(false);
        _best = _ener.print(_ener.N_blk(),0,false);
    };
    double print(int blk,int prop,bool prt)override;
    vec increase(vec pos)override;
    void adjust(int prop)override;
    
private:
    vec _min;
    double _best;
    energy& _ener;
    metropolis& _metro;
    double _t, _div;
};


#endif /* datablocking_h */
