//
//  metropolis.hpp
//
//
//  Created by Giovanni Galafassi on 08/08/24.
//
#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdio.h>

#include "population.h"

using namespace std;
using namespace arma;

void gene::Naming(string n){
    _name = n;
};
void gene::SetId(int i){
    _label = i;
};
void gene::Placing(vec x){
    _pos = x;
};
void gene::Init(string name,int i,vec x){
    _name = name;
    _label = i;
    _pos = x;
};
int gene::Id(){return _label;};
vec gene::Pos()const{return _pos;};


gene& element::GetG(int i){
    return _dna[i];
};
int element::N(){
    return _dim;
};
bool element::Check(){
    int count = 0;
    for(int i=0;i<_dim;i++){
        for(int j=i+1;j<_dim;j++){
            if(_dna[i].Id() == _dna[j].Id()) return false;
        }
        count += _dna[i].Id();
    }
    if (count == _sum) return true;
    else return false;
};
double element::L2()const{
    double l2 = 0.;
    for(int i=0;i<_dim-2;i++){
        l2 += norm(_dna[i].Pos() - _dna[i+1].Pos());
    }
    l2 += norm(_dna[_dim -1].Pos() - _dna[0].Pos());
    return l2;
};
void element::SetL2(){
    _l2 = L2();
};
int element::PBC(int i,int dim){
    return i%dim;
};

vector <gene> element::Cut(int i, int j){
    if(j == -1) j = _dim;
    vector <gene> app(j - i);
    for(int k=i;k<j;k++) app[k-i] = _dna[k];
    return app;
};
vector <gene> element::Cut(vector <gene> a,int i, int j){
    if(j == -1) j = _dim;
    vector <gene> app(j - i);
    for(int k=i;k<j;k++) app[k-i] = a[PBC(k,a.size())];
    return app;
};
void element::Rebuild(vector <gene> r,int i){
    for(int k=i;k<_dim;k++) _dna[k] = r[k-i];
};
void element::Shift(int n, int m, int i){
    i--;
    vector <gene> rap = Cut(1);
    vector <gene> uno = Cut(rap,i,i+m);
    vector <gene> due = Cut(rap,i+m,i+m+n);
    
    for(int k=i;k<i+n;k++) rap[PBC(k,_dim-1)] = uno[k];
    for(int k=i+n;k<i+n+m;k++) rap[PBC(k,_dim-1)] = due[k];
    //Rebuild(rap);
    for(int k=1;k<_dim;k++) _dna[k] = rap[k-1];
};

void element::RndSwap(){
    int i,j;
    do{
        i = int(_rnd.Rannyu(1,_dim));
        j = int(_rnd.Rannyu(1,_dim));
    }while( i == j );
    gene appo = _dna[i];
    _dna[i] = _dna[j];
    _dna[j] = appo;
};
void element::RndSwap(int m){
    int i,j;
    do{
        i = int(_rnd.Rannyu(1,_dim));
        j = int(_rnd.Rannyu(1,_dim));
    }while( abs(i-j) < m );
    i--; j--;
    vector <gene> rap = Cut(1);
    vector <gene> uno = Cut(rap,i,i+m);
    vector <gene> due = Cut(rap,j,j+m);
    
    for(int k=i;k<i+m;k++) rap[PBC(k,_dim-1)] = due[k];
    for(int k=j;k<j+m;k++) rap[PBC(k,_dim-1)] = uno[k];
    Rebuild(rap);
};
void element::Inv(int i, int j){
    i--;j--;
    vector <gene> rap = Cut(1);
    vector <gene> uno = Cut(rap,i,j);
    
    for(int k=i;k<j;k++) rap[PBC(k,_dim-1)] = uno[j-k-1];
    Rebuild(rap);
};
void element::RndInv(){
    int i,j;
    do{
        i = int(_rnd.Rannyu(1,_dim));
        j = int(_rnd.Rannyu(1,_dim));
    }while( i >= j );
    Inv(i,j);
};

element& population::GetEl(int i){
    return _pop[i];
};
element& population::RanEl(){
    return _pop[int(_rnd.Rannyu() * _n)];
};
void population::Circle(int length){
    cout<<"Preparing Adamo"<<endl;
    element adamo(length,_rnd);
    for(int i=0;i<length;i++){
        double b = _rnd.Rannyu() * 2. * M_PI;
        vec pos = {sin(b),cos(b)};
        adamo(i).Placing(pos);
        adamo(i).SetId(i+1);
        Progress_Bar(i,length);
    }
    cout<<endl<<"Adamo created"<<endl;
    for(int j=0;j<_n;j++) _pop.push_back(adamo);
};
void population::Spread(){
    cout<<"Spreading population"<<endl;
    for(int j=0;j<_n;j++){
        for(int k=0;k<5;k++){
            //_pop[j].RndSwap();
            //_pop[j].RndSwap(2);
            _pop[j].Shift(10,3,30);
        }
        Progress_Bar(j,_n);
    }
    cout<<endl<<"Done"<<endl;
};
int population::Check(){
    int ok = 0;
    for(int i=0;i<_n;i++) if( ! _pop[i].Check() ) ok = i+1;
    return ok;
};
void population::Xover(int i, int j){
    vector <gene> uno = _pop[i].Cut(1);
    vector <gene> due = _pop[j].Cut(1);
    int k = int(_rnd.Rannyu() * (_pop[i].N() - 1));
    
    vector <gene> c1 = _pop[i].Cut(uno,k);
    vector <gene> c2 = _pop[j].Cut(due,k);
    
    Sort(c1,due);
    Sort(c2,uno);
    
    _pop[i].Rebuild(c1,k);
    _pop[j].Rebuild(c2,k);
};
int population::Search(gene x, vector <gene> dad){
    int index = -1;
    for(int i=0;i<dad.size();i++){
        if(x.Id() == dad[i].Id()) index = i;
    }
    return index;
};
void population::Sort(vector <gene> guy, vector <gene> dad){
    sort(guy.begin(),guy.end(),[this,&dad](const gene& a, const gene& b){
        return this->Search(a,dad) < this->Search(b,dad);
    });
};
void population::Sort(){
    sort(_pop.begin(),_pop.end(),[](const element& a, const element& b){
        return a.L2() < b.L2();
    });
};
void population::Print(int i){
    ofstream fout("element_"+to_string(i)+".txt");
    for(int j=0;j<_pop[i].N();j++){
        fout<<_pop[i](j).Id();
        for(int k=0;k<_pop[i](j).Pos().n_elem;k++){
            fout<<"\t"<<_pop[i](j).Pos()(k);
        }
        fout<<endl;
    }
    fout.close();
};
