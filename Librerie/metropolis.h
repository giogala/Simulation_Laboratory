//
//  metropolis.hpp
//
//
//  Created by Giovanni Galafassi on 09/05/24.
//  Modified on 16/05/24.
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
        int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
        ifstream Primes(rnd_dir+"/Primes");
        Primes >> p1 >> p2 ;
        Primes.close();
        int seed[4]; // Read the seed of the RNG
        ifstream Seed(rnd_dir+"/seed.in");
        Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        _rnd.SetRandom(seed,p1,p2);
        
        //acceptance things
        tried=0;
        accepted=0;
        sigma=distr.GetParams(0);
        _rnd.Gauss(0.,sigma);
        _rnd.Gauss(0.,sigma);
    };
    ~metropolis(){;};
    
    bool metro(double prob){
        bool accept = false;
        if(_rnd.Rannyu() < prob) accept = true;
        return accept;
    };
    
    vec move(vec init){
        tried++;
        vec fin = init;
        if(_type=="gauss"){
            for(int i=0;i<fin.n_elem;i++) fin(i) = _rnd.Gauss(0.,sigma);
            //fin = {_rnd.Gauss(0.,sigma),_rnd.Gauss(0.,sigma),_rnd.Gauss(0.,sigma)};
            fin = fin + init;
        }
        else if(_type=="unif"){
            for(int i=0;i<fin.n_elem;i++) fin(i) = _rnd.Rannyu(-sigma,sigma);
            //fin = {_rnd.Rannyu(-sigma,sigma),_rnd.Rannyu(-sigma,sigma),_rnd.Rannyu(-sigma,sigma)};
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
    
    double trans_prob(vec init, vec fin){
        double t1 = distr.eval(fin);
        double t2 = distr.eval(init);
        return min(1.,t1/t2);
    };

    
    void adjust(){
        // riscalo la \sigma della funzione di transizione in modo da avere accettanza al 50%
        sigma*=((double)accepted)/((double)tried*0.5);
        tried=0; accepted=0;
    }
    
    funzione& GetDistr(){
        return distr;
    }
    
protected:
    funzione& distr;
    string _type;
    double sigma;
    Random _rnd;
    int tried;
    int accepted;
};
#endif /* metropolis.h */

