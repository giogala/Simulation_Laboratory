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
    funzione():par(1){;};
    virtual ~funzione(){;};
    
    virtual double eval(double a) const =0;
    virtual double eval(vec a) const =0;
    virtual double der(vec a) const =0;
    virtual string whoamI() =0;
    
    vec radial(vec a)const{
        double r = sqrt(dot(a,a));
        double rp = pow(a(0),2)+pow(a(1),2);
        double phi = 2*atan(a(1)/(sqrt(rp)+a(0)));
        double theta = acos(a(2)/r);
        vec b={r,theta,phi};
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
    
    void SetParams(vec a){
        if(a.n_elem == n_par) par = a;
        else{
            for(int i=0;i<(int)min((double)a.n_elem,(double)n_par);i++) par(i)=a(i);
        }
    };
    
    void SetParam(int i, double a){
        if(i<n_par) par(i)=a;
        else{
            cerr<<"Error: can't reach "<<i<<"-th parameter: only "<<n_par<<" are given"<<endl;
            exit(EXIT_FAILURE);
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
        par(0)=1.;//0.0529; //a_0 raggio di Bohr in [nm]
        rad=false;
    };
    ~hydro_100(){;};
    
    double eval(double a)const{
        double r = a;
        return pow(par(0),-3)*exp(-2*r/par(0))/M_PI;
    };
    double eval(vec a)const{
        double r = sqrt(dot(a,a));
        return pow(par(0),-3)*exp(-2*r/par(0))/M_PI;
    };
    double der(vec a) const{return 0;};
    string whoamI(){
        return "hydro_100";
    };
    
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
    
    double eval(double a)const{
        return -1;
    }
    
    double eval(vec a)const{
        vec r = a;
        if(rad) r = radial(r);
        return pow(par(0),-5)*r(0)*r(0)*exp(-r(0)/par(0))*pow(cos(r(1)),2)/(32*M_PI);
    };
    double der(vec a) const{return 0;};
    string whoamI(){
        return "hydro_210";
    };
};

class hat: public funzione{
public:
    hat(double mu,double sigma){
        n_par=2;
        par.set_size(n_par);
        par = {sigma,mu};
        rad=false;
    };
    ~hat(){;};
    
    double eval(vec a)const{
        if(a.n_elem == 1){
            double nor = 1./(pow(M_PI,1./4.)*pow(par(0)*(1.+2.*exp(-pow(par(0)/par(1),2))),1./2.));
            double uno = exp(-pow((a(0)-par(1))/par(0),2));
            double due = exp(-pow((a(0)+par(1))/par(0),2));
            double tre = 2.*exp(-(a(0)*a(0)+par(1)*par(1))/(par(0)*par(0)));
            return nor*(uno+due+tre);
        }
        else{
            cerr<<"Error: only 1D vectors supported for hat function"<<endl;
            exit(EXIT_FAILURE);
            return -1;
        }
    };
    double eval(double r)const{
        double nor = 1./(pow(M_PI,1./4.)*pow(par(0)*(1.+2.*exp(-pow(par(0)/par(1),2))),1./2.));
        double uno = exp(-pow((r-par(1))/par(0),2));
        double due = exp(-pow((r+par(1))/par(0),2));
        double tre = 2.*exp(-(r*r+par(1)*par(1))/(par(0)*par(0)));
        return nor*(uno+due+tre);
    };
    double der(vec a)const{
        if(a.n_elem == 1){
            double x = a(0);
            double mu = par(1);
            double s = par(0)*par(0);
            double num = pow((x-mu)/s,2)*exp(-pow(x-mu,2)/(2*s)) + pow((x+mu)/s,2)*exp(-pow(x+mu,2)/(2*s));
            double den = exp(-pow(x-mu,2)/(2*s)) + exp(-pow(x+mu,2)/(2*s));
            double pot = pow(x,4) - 2.5*pow(x,2);
            return 0.5*(1./s - num/den) + pot;
        }
        else{
            cerr<<"Error: only 1D vectors supported for hat function"<<endl;
            exit(EXIT_FAILURE);
            return -1;
        }
    };
    string whoamI(){
        return "hat";
    }
};

class blotzmann: public funzione{
public:
    boltzmann(double temp){
        n_par=1;
        par.set_size(n_par);
        par = {temp};
        rad=false;
    };
    ~boltzmann(){;};
    
    double eval(vec a)const{
        if(a.n_elem == 1){
            double e = a(0);
            return exp(-e/par(0));
        }
        else{
            cerr<<"Error: only 1D vectors supported for Boltzmann distribution"<<endl;
            exit(EXIT_FAILURE);
            return -1;
        }
    };
    double eval(double r)const{
        return exp(-r/par(0));
    };
    double der(vec a)const{return 0;};
    string whoamI(){
        return "boltzmann";
    }
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
    
    double eval(vec a)const{
        vec r = a;
        if(rad) r = radial(r);
        vec c = {par(1),par(2),par(3)};
        //cerr<<c<<endl;
        return pow(2*M_PI*par(0)*par(0),-(double)(3/2))*exp(-dot(r-c,r-c)/(2*par(0)*par(0)));
    };
    double eval(double r)const{
        return exp(-r*r/(2*par(0)*par(0)))/sqrt(pow(2*M_PI*par(0)*par(0),3));
    };
    string whoamI(){
        return "gauss";
    }
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
    
    double eval(vec a)const{return 1.;};
    double eval(double a)const{return 1.;};
    string whoamI(){
        return "unif";
    };
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
