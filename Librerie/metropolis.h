//
//  metropolis.hpp
//
//
//  Created by Giovanni Galafassi on 16/05/24.
//  Modified on 09/05/24.
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

class metropolis:{
public:
    metropolis(){;};
    metropolis(funzione p,string t, vec pos){
        distr=p;
        type=t;
        if(type=='gauss'){
            gauss3d tra();
            trans = tra;
        }
        if(type=='unif'){
            unif tra();
            trans = tra;
        }
        else{
            cerr<<"Erroe: Unknown transition distribution type given"<<endl;
            exit(EXIT_FAILURE);
        }
        pos_0=pos;
        int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
        ifstream Primes("Primes");
        Primes >> p1 >> p2 ;
        Primes.close();
        int seed[4]; // Read the seed of the RNG
        ifstream Seed("seed.in");
        Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        _rnd.SetRandom(seed,p1,p2);
        tried=0;
        accepted=0;
        sigma=1.;
    };
    ~metropolis(){;};
    
    bool metro(double prob){
        bool accept = false;
        if(_rnd.Rannyu() < prob) accept = true;
        return accept;
    };
    
    vec move(vec init){
        tried++; accepted++;
        if(type=='gauss') vec fin = init + {_rnd.Gauss(0,sigma),_rnd.Gauss(0,sigma),_rnd.Gauss(0,sigma)};
        if(type=='unif') vec fin = init + {_rnd.Rannyu(-sigma,sigma),_rnd.Rannyu(-sigma,sigma),_rnd.Rannyu(-sigma,sigma)};
        
        if(! metro(trans_prob(distr,init,fin))) {
            fin = init;
            accepted--;
        }
        return init;
    };
    
    double trans_prob(vec init, vec fin){
        trans.SetParams(fin); //centro la funzione di transizione in pos_f
        double t1 = trans.eval(init)*distr.eval(fin);
        trans.SetParams(init); //centro la funzione di transizione in pos_i
        double t2 = trans.eval(fin)*distr.eval(init);
        return min(1,t1/t2);
    };
    // da mettere a posto la funzione di transizione. Devo decidere se porre il centro della distribuzione come variabile della funzione eval o come data_membro della classe. Sarei forse più per la seconda ma i parametri di "funzione" sono implementati come un vec e se il mio parametro 0 fosse il centro avrei che vec(0) sarebbe a sua volta un 3d-vec mentre vec(1), ovvero la \sigma, sarebbe un double. Non penso che vec possa contenere diversi tipi di dati quindi come alternativa potrei avere una vec siffatto: (x_c,y_c,z_c,\sigma)
    
    // ora, questa cosa, dovrebbe essere a posto
    
    void adjust(){
        // riscalo la \sigma della funzione di transizione in modo da avere accettanza al 50%
        if(type=='gauss') simga=GetParam(3)*0.5*(double)(tried/accepted)
        if(type=='unif') simga=GetParam(0)*0.5*(double)(tried/accepted)
        trans.SetParam(i,sigma);
        tried=0; accepted=0;
    }
    
protected:
    funzione distr;
    funzione trans;
    string type;
    double sigma;
    vec pos_0;
    Random _rnd;
    int tried;
    int accepted;
};
#endif /* metropolis.h */
/*hydro(int a,int b, int c){
    par(0)=a; par(1)=b; vec(2)=c;
};
coseno(const coseno& s){
    m_a=s.geta(); m_b=s.getb(); m_c=s.getc();
};*/
/*void SetPar(int i, int p){
    par(i)=p;
};

int GetPar(int i)const{
    return par(i);
};

//double xv();
//double yv();

//coseno& operator=(const coseno&);
//coseno& operator=(const seno&);
*/
