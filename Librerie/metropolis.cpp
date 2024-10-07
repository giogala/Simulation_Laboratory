//
//  metropolis.hpp
//
//
//  Created by Giovanni Galafassi on 09/05/24.
//
#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdio.h>

#include "metropolis.h"

using namespace std;
using namespace arma;


bool metropolis::metro(double prob){
    bool accept = false;
    if(_rnd.Rannyu() < prob) accept = true;
    return accept;
};

vec metropolis::move(vec init){
    tried++;
    vec fin = init;
    if(_type=="gauss"){
        for(int i=0;i<fin.n_elem;i++) fin(i) = _rnd.Gauss(0.,sigma);
        fin = fin + init;
    }
    else if(_type=="unif"){
        for(int i=0;i<fin.n_elem;i++) fin(i) = _rnd.Rannyu(-sigma,sigma);
        fin = fin + init;
    }
    else{
        cerr<<"Error: Unknown transition distribution type given"<<endl;
        exit(EXIT_FAILURE);
    }
    if(metro(trans_prob(init,fin))) accepted++;
    else fin = init;
    return fin;
};

double metropolis::trans_prob(vec init, vec fin){
    double t1 = distr.eval(fin);
    double t2 = distr.eval(init);
    return min(1.,t1/t2);
};

void metropolis::adjust(){
    // riscalo la \sigma della funzione di transizione in modo da avere accettanza al 50%
    sigma*=((double)accepted)/((double)tried*0.5);
    tried=0; accepted=0;
};

void metropolis::SetS(double s){
    sigma = s;
};

funzione& metropolis::GetDistr(){
    return distr;
};
