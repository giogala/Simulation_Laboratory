//
//  funzione.cpp
//
//
//  Created by Giovanni Galafassi on 10/11/21.
//
#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdio.h>
#include "funzione.h"

using namespace std;
using namespace arma;


vec funzione::radial(vec a)const{
    double r = sqrt(dot(a,a));
    double rp = pow(a(0),2)+pow(a(1),2);
    double phi = 2*atan(a(1)/(sqrt(rp)+a(0)));
    double theta = acos(a(2)/r);
    vec b={r,theta,phi};
    return b;
}
double funzione::GetParams(int i){
    if(i<n_par) return par(i);
    else{
        cerr<<"Error: can't reach "<<i<<"-th parameter: only "<<n_par<<" aregiven"<<endl;
        exit(EXIT_FAILURE);
        return -1;
    }
};

void funzione::SetRadial(bool tof){
    rad=tof;
};

void funzione::SetParams(vec a){
    if(a.n_elem == n_par) par = a;
    else{
        for(int i=0;i<(int)min((double)a.n_elem,(double)n_par);i++) par(i)=a(i);
    }
};

void funzione::SetParam(int i, double a){
    if(i<n_par) par(i)=a;
    else{
        cerr<<"Error: can't reach "<<i<<"-th parameter: only "<<n_par<<" aregiven"<<endl;
        exit(EXIT_FAILURE);
    }
};



double hydro_100::eval(double a)const{
    double r = a;
    return pow(par(0),-3)*exp(-2*r/par(0))/M_PI;
};
double hydro_100::eval(vec a){
    double r = sqrt(dot(a,a));
    return pow(par(0),-3)*exp(-2*r/par(0))/M_PI;
};
double hydro_100::der(vec a) const{return 0;};
string hydro_100::whoamI(){
    return "hydro_100";
};



double hydro_210::eval(double a)const{
    return -1;
};

double hydro_210::eval(vec a){
    vec r = a;
    if(rad) r = radial(r);
    return pow(par(0),-5)*r(0)*r(0)*exp(-r(0)/par(0))*pow(cos(r(1)),2)/(32*M_PI);
};
double hydro_210::der(vec a) const{return 0;};
string hydro_210::whoamI(){
    return "hydro_210";
};



double hat::eval(vec a){
    if(a.n_elem == 1){
        //double nor = 1./(pow(M_PI,1./4.)*pow(par(0)*(1.+2.*exp(-pow(par(0)/par(1),2))),1./2.));
        double x = a(0);
        double mu = par(1);
        double s = par(0);
        double uno = exp( -pow(x - mu,2)/(s*s) );
        double due = exp( -pow(x + mu,2)/(s*s) );
        double tre = 2. * exp( -pow(x - mu,2)/(2.*s*s) ) * exp( -pow(x + mu,2)/(2.*s*s) );
        return uno+due+tre;
    }
    else{
        cerr<<"Error: only 1D vectors supported for hat function"<<endl;
        exit(EXIT_FAILURE);
        return -1;
    }
};
double hat::eval(double r)const{
    double nor =1./(pow(M_PI,1./4.)*pow(par(0)*(1.+2.*exp(-pow(par(0)/par(1),2))),1./2.));
    double uno = exp(-pow((r-par(1))/par(0),2));
    double due = exp(-pow((r+par(1))/par(0),2));
    double tre = 2.*exp(-(r*r+par(1)*par(1))/(par(0)*par(0)));
    return nor*(uno+due+tre);
};
double hat::der(vec a)const{
    if(a.n_elem == 1){
        long double x = a(0);
        long double mu = par(1);
        long double s = par(0);
        
        double kin = -0.5 * ((pow(x - mu,2)/pow(s,4) - 1./pow(s,2)) *exp(-(pow(x - mu,2)/(2.*s*s))) + (pow(x + mu,2)/pow(s,4) - 1./pow(s,2)) *exp(-(pow(x + mu,2)/(2.*s*s))))/(exp(-(pow(x - mu,2)/(2.*s*s))) + exp(-(pow(x + mu,2)/(2.*s*s))));
        double pot = pow(x, 4) - 2.5 * pow(x, 2);
        return kin + pot;
    }
    else{
        cerr<<"Error: only 1D vectors supported for hat-function"<<endl;
        exit(EXIT_FAILURE);
        return -1;
    }
};
string hat::whoamI(){
    return "hat";
};



double gauss3d::eval(vec a){
    vec r = a;
    if(rad) r = radial(r);
    vec c = {par(1),par(2),par(3)};
    return pow(2*M_PI*par(0)*par(0),-(double)(3/2))*exp(-dot(r-c,r-c)/(2*par(0)*par(0)));
};
double gauss3d::eval(double r)const{
    return exp(-r*r/(2*par(0)*par(0)))/sqrt(pow(2*M_PI*par(0)*par(0),3));
};
string gauss3d::whoamI(){
    return "gauss";
};



double unif::eval(vec a){return 1.;};
double unif::eval(double a)const{return 1.;};
string unif::whoamI(){ return "unif";};

