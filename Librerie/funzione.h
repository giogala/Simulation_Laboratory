//
//  funzione.hpp
//  
//
//  Created by Giovanni Galafassi on 10/11/21.
//
#define _USE_MATH_DEFINES
#ifndef __funzione_h_
#define __funzione_h_
#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdio.h>
//#include "datablocking.h"


using namespace std;
using namespace arma;

class funzione{
public:
    funzione():par(1){;};
    virtual ~funzione(){;};
    
    virtual double eval(double a) const =0;
    virtual double eval(vec a) =0;
    virtual double der(vec a) const =0;
    virtual string whoamI() =0;
    
    double operator()(vec x) {return eval(x);}
    vec radial(vec a)const;
    double GetParams(int i);
    void SetRadial(bool tof);
    void SetParams(vec a);
    void SetParam(int i, double a);
    
protected:
    int n_par;
    vec par;
    bool rad;
};



class hydro_100: public funzione{
public:
    hydro_100(){
        n_par=1;
        par.set_size(n_par);
        par(0)=1.;//0.0529; //a_0 raggio di Bohr in [nm]
        rad=false;
    };
    ~hydro_100(){;};
    
    double eval(double a)const;
    double eval(vec a);
    double der(vec a) const;
    string whoamI();
};

class hydro_210: public funzione{
public:
    hydro_210(){
        n_par=1;
        par.set_size(n_par);
        par(0)=1.;//0.0529; //a_0 raggio di Bohr in [nm]
        rad=true;
    };
    ~hydro_210(){;};
    
    double eval(double a)const;
    double eval(vec a);
    double der(vec a) const;
    string whoamI();
};

class hat: public funzione{
public:
    hat(double sigma,double mu){
        n_par=2;
        par.set_size(n_par);
        par = {sigma,mu};
        rad=false;
    };
    ~hat(){;};
    
    double eval(vec a);
    double eval(double r)const;
    double der(vec a)const;
    string whoamI();
};

class energy;

class boltzmann: public funzione{
public:
    boltzmann(double temp,energy& f):_f(f){
        n_par=1;
        par.set_size(n_par);
        par = {temp};
        rad=false;
    };
    ~boltzmann(){;};
    
    double eval(vec a);
    double eval(double r)const;
    double der(vec a)const;
    string whoamI();
    
private:
    energy& _f;
    double _en;
};

class gauss3d: public funzione{
public:
    gauss3d(){
        n_par=4;
        par.set_size(n_par);
        par = {1.,0.,0.,0.};
        rad=false;
    };
    ~gauss3d(){;};
    
    double eval(vec a);
    double eval(double r)const;
    string whoamI();
};

class unif: public funzione{
public:
    unif(){
        n_par=4;
        par.set_size(n_par);
        par = {1.,0.,0.,0.};
        rad=false;
    };
    ~unif(){;};
    
    double eval(vec a);
    double eval(double a)const;
    string whoamI();
};


#endif /* funzione.h */

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
