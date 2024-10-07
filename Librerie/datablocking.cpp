//
//  datablocking.h
//
//
//  Created by Giovanni Galafassi on 02/10/24.
//
#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdio.h>
#include "datablocking.h"

using namespace std;
using namespace arma;


void datablocking::blocks(bool prt){
    for(int i=0;i<_nblock;i++){
        blk_av.zeros();
        _state = increase(_state);
        if (prt) Progress_Bar(i+1,_nblock);
        for(int j=0;j<_nprop;j++) {
            averages(j);
            print(i+1,j,prt);
        }
    }
};

void datablocking::averages(int prop){
    _av(prop) = blk_av(prop)/double(_nstep);
    adjust(prop);
    glb_av(prop) += _av(prop);
    glb_2(prop) += _av(prop) * _av(prop);
};

double datablocking::print(int blk,int prop,bool prt){
    average  = _av(prop);
    sum_average = glb_av(prop);
    sum_ave2 = glb_2(prop);
    if(prt){
        ofstream fout;
        fout.open(_f[prop]+".txt",ios::app);
        fout << blk
        << "\t" << average
        << "\t" << sum_average/double(blk)
        << "\t" << this->error(sum_average, sum_ave2, blk) << endl;
        fout.close();
    }
    return sum_average/double(blk);
};

double datablocking::error(double acc, double acc2, int blk){
    if(blk <= 1) return 0.0;
    else return sqrt( fabs(acc2/double(blk) - pow( acc/double(blk) ,2) )/double(blk) );
};
void datablocking::reset(){
    _av.zeros();
    blk_av.zeros();
    glb_av.zeros();
    glb_2.zeros();
};
int datablocking::N_blk(){return _nblock;};

vec datablocking::pos(){return _state;};

vec prova::increase(vec pos){
    for(int k=0;k<_nstep;k++){
        pos = _metro.move(pos);
        blk_av(0) += sqrt(dot(pos,pos));
    }
    return pos;
};
void prova::adjust(int prop){
    if(prop==0)_metro.adjust();
};


vec energy::increase(vec pos){
    for(int k=0;k<_nstep;k++){
        pos = _metro.move(pos);
        blk_av(0) += _metro.GetDistr().der(pos);
        if(_prt){
            ofstream fout;
            fout.open("position.txt",ios::app);
            fout<<pos(0)<<endl;
            fout.close();
        }
    }
    return pos;
};
void energy::adjust(int prop){
    if(prop==0)_metro.adjust();
};
metropolis& energy::GetMetro(){return _metro;};

double annealing::print(int blk,int prop,bool prt){
    average  = _av(prop);
    sum_average = glb_av(prop);
    sum_ave2 = glb_2(prop);
    if(prt){
        ofstream fout;
        fout.open(_f[prop]+".txt",ios::app);
        if(prop == 0){
            fout << _t * _div
            << "\t" << _min(0)
            << "\t" << _min(1)
            << "\t" << _best << endl;
        } else {
            fout << blk
            << "\t" << average
            << "\t" << sum_average/double(blk)
            << "\t" << this->error(sum_average, sum_ave2, blk) << endl;
        }
        fout.close();
    }
    return sum_average/double(blk);
};

vec annealing::increase(vec pos){
    for(int k=0;k<_nstep;k++){
        pos = _metro.move(pos);
        double ef = _metro.GetDistr().der(pos);
        if( ef <= _best){
            _best = ef;
            _min = pos;
        }
        blk_av(1) += ef;
    }
    blk_av(2) = pos(0) * _nstep;
    blk_av(3) = pos(1) * _nstep;
    return pos;
};

void annealing::adjust(int prop){
    if(prop==0){
        _t /= _div;
        _metro.GetDistr().SetParam(0,_t);
        //_metro.adjust();
    }
};

double boltzmann::eval(vec a){
    _f.reset();
    _f.GetMetro().GetDistr().SetParams(a);
    _f.blocks(false);                       // eseguo data_blocking per \mu e \sigma fissati
    int blk = _f.N_blk();
    double e = _f.print(blk,0,false);       // estraggo il valore finale di E
    _en = e;
    _f.reset();                             // resetto il data_blocking
    return exp(-e/par(0));                   // calcolo e(a) l'energia in funzione dalla "posizione"
    
};
double boltzmann::eval(double r)const{
    return exp(-r/par(0));
};
double boltzmann::der(vec a)const{return _en;};
string boltzmann::whoamI(){
    return "boltzmann";
};
