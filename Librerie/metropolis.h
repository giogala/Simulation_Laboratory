//
//  metropolis.hpp
//
//
//  Created by Giovanni Galafassi on 09/05/24.
//  
#define _USE_MATH_DEFINES
#ifndef __metropolis_h_
#define __metropolis_h_
#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdio.h>
#include "funzione.h"
#include "random.h"

using namespace std;
using namespace arma;

class metropolis{
public:

    metropolis(funzione& p,string t, string rnd_dir):distr(p){
        // transition distribution
        _type=t;
        
        // Random generator
        _rnd.initRnd(rnd_dir);
        //acceptance things
        tried=0;
        accepted=0;
        sigma=distr.GetParams(0);
        _rnd.Gauss(0.,sigma);
        _rnd.Gauss(0.,sigma);
    };
    ~metropolis(){;};
    
    bool metro(double prob);
    vec move(vec init);
    double trans_prob(vec init, vec fin);
    void adjust();
    void SetS(double s);
    funzione& GetDistr();
    
protected:
    funzione& distr;
    string _type;
    double sigma;
    Random _rnd;
    int tried;
    int accepted;
};
#endif /* metropolis.h */

