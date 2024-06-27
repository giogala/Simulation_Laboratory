//
//  funzione.hpp
//  
//
//  Created by Giovanni Galafassi on 10/11/21.
//  Modified on 09/05/24.
#define _USE_MATH_DEFINES
#ifndef __funzione_h_
#define __funzione_h_
#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdio.h>

using namespace std;
using namespace arma;

class funzione{
public:
    funzione(){;};
    virtual ~funzione(){;};
    
    virtual double eval(double a) const =0;
    virtual double eval(vec a) const =0;
    
    
    vec radial(vec a){
        double r = sqrt(a*a);
        double rp = pow(a(0),2)+pow(a(1),2);
        double phi = 2*atan(a(1)/(sqrt(rp)+a(0)));
        double theta = 2*atan(rp/(sqrt(rp*rp+a(2)*a(2))+a(2)));
        vec b=(r,theta,phi);
        return b;
    }
    double GetParams(int i){
        if(i<n_par) return par(i);
        else{
            cerr<<"Error: can't reach "<<i<<"-th parameter: only "<<n_par<<" are given"<<endl;
            exit(EXIT_FAILURE);
            return -1;
        }
    };
    
    void SetRadial(bool tof){
        rad=tof;
    };
    
    void SetParams(vec a)const{
        if(a.n_elem == n_par) par = a;
        else{
            for(int i=0;i<min(a.n_elem,n_par);i++) n_par(i)=a(i);
        }
    };
    
    void SetParam(int i, double a)const{
        if(i<n_par) par(i)=a;
        else{
            cerr<<"Error: can't reach "<<i<<"-th parameter: only "<<n_par<<" are given"<<endl;
            exit(EXIT_FAILURE);
            return -2;
        }
    };
    
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
        par(0)=0.0529 //a_0 raggio di Bohr in [nm]
    };
    ~hydro_100(){;};
    
    double eval(double a)const{
        double r = a;
        return pow(par(0),-3./2.)*exp(-r/par(0))/sqrt(M_PI);
    };
    double eval(vec a)const{
        double r = sqrt(a*a);
        return pow(par(0),-3./2.)*exp(-r/par(0))/sqrt(M_PI);
    };
    
};

class hydro_210: public funzione{
public:
    hydro_210(){
        n_par=1;
        par.set_size(n_par);
        par(0)=0.0529 //a_0 raggio di Bohr in [nm]
    };
    ~hydro_210(){;};
    
    double eval(vec a)const{
        if(rad) vec r = radial(a);
        else: vec r = a;
        return pow(par(0),-5./2.)*sqrt(2)*r(0)*exp(-r(0)/(2*par(0)))*cos(r(1))/(8*sqrt(M_PI));
    };
};

class gauss3d: public funzione{
public:
    gauss3d(){
        n_par=4;
        par.set_size(n_par);
        par = {0.,0.,0.,1.};
    };
    ~gauss3d(){;};
    
    double eval(vec a)const{
        if(rad) vec r = radial(a);
        else: vec r = a;
        vec c = {par(0),par(1),par(2)};
        return pow(2*M_PI*par(3)*par(3),-3/2)*exp(-(r-c)*(r-c)/(2*par(3)*par(3)));
    };
};

class unif: public funzione{
    unif(){
        n_par=1;
        par.set_size(n_par);
        par = {1.};
    };
    ~unif(){;};
    
    double eval(vec a)const{
        return 1.;
    };
}


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
