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
int gene::Search(vector <gene> dad){
    int index = -1;
    for(int i=0;i<dad.size();i++){
        if(_label == dad[i].Id()) index = i;
    }
    return index;
};


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
    return (i % dim + dim) % dim;
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
void element::Rebuild(vector <gene> zip,vector <gene> r,int i){
    for(int k=i;k<zip.size();k++) zip[k] = r[k-i];
};
void element::Shift(int n, int m, int i){
    i--;
    vector <gene> rap = Cut(1);
    vector <gene> uno = Cut(rap,i,i+m);
    vector <gene> due = Cut(rap,i+m,i+m+n);
    
    for(int k=i;k<i+n;k++) rap[PBC(k,_dim-1)] = due[k-i];
    for(int k=i+n;k<i+n+m;k++) rap[PBC(k,_dim-1)] = uno[k-i-n];
    Rebuild(rap);
};
void element::Swap(int i,int j,int m){
    i--; j--;
    vector <gene> rap = Cut(1);
    vector <gene> uno = Cut(rap,i,i+m);
    vector <gene> due = Cut(rap,j,j+m);
    
    for(int k=i;k<i+m;k++) rap[PBC(k,_dim-1)] = due[k-i];
    for(int k=j;k<j+m;k++) rap[PBC(k,_dim-1)] = uno[k-j];
    Rebuild(rap);
};
void element::Inv(int i, int j){
    i--;j--;
    vector <gene> rap = Cut(1);
    vector <gene> uno = Cut(rap,i,j);
    
    for(int k=i;k<j;k++) rap[PBC(k,_dim-1)] = uno[j-k-1];
    Rebuild(rap);
};

element& population::GetEl(int i){
    return _pop[i];
};
element& population::RanEl(){
    return _pop[int(pow(_rnd.Rannyu(),5) * _n)];
};
int population::Select(){
    return int(pow(_rnd.Rannyu(),5) * _n);
};
void population::Circle(int length){
    element adamo(length);
    for(int i=0;i<length;i++){
        double b = _rnd.Rannyu() * 2. * M_PI;
        vec pos = {sin(b),cos(b)};
        adamo(i).Placing(pos);
        adamo(i).SetId(i+1);
    }
    for(int j=0;j<_n;j++) _pop.push_back(adamo);
    Spread();
};
void population::Square(int length){
    element adamo(length);
    for(int i=0;i<length;i++){
        double x = _rnd.Rannyu();
        double y = _rnd.Rannyu();
        vec pos = {x,y};
        adamo(i).Placing(pos);
        adamo(i).SetId(i+1);
    }
    for(int j=0;j<_n;j++) _pop.push_back(adamo);
    Spread();
};
void population::Spread(){
    for(int j=1;j<_n;j++){
        for(int k=0;k<3;k++){
            RndSwap(j);
            RndSwap(j,2);
            RndSwap(j,4);
            RndInv(j);
        }
    }
};
void population::Check(){
    for(int i=0;i<_n;i++) if( ! _pop[i].Check() ) cout<<i<<endl;
};
void population::RndSwap(int k, int m){
    int dim = _pop[k].N();
    int i,j;
    do{
        i = int(_rnd.Rannyu(1,dim));
        j = int(_rnd.Rannyu(1,dim));
    }while( _pop[k].PBC(i,dim-1) - _pop[k].PBC(j+m,dim-1) or _pop[k].PBC(i,dim-1) > _pop[k].PBC(j-m,dim-1));
    _pop[k].Swap(i,j,m);
};
void population::RndInv(int k){
    int dim = _pop[k].N();
    int i,j;
    do{
        i = int(_rnd.Rannyu(1,dim));
        j = int(_rnd.Rannyu(1,dim));
    }while( i >= j );
    _pop[k].Inv(i,j);
};
void population::RndShift(int k){
    int m = _rnd.Rannyu(1, _pop[k].N()-1);
    int n = _rnd.Rannyu(1, _pop[k].N()-m);
    int i = _rnd.Rannyu(1,_pop[k].N());
    _pop[k].Shift(n,m,i);
};
void population::Xover(int i, int j){
    int l,m;
    do{
        l = Select();
        m = Select();
    } while( l == m);
    
    vector <gene> uno = _pop[l].Cut(1);
    vector <gene> due = _pop[m].Cut(1);
    int d = _pop[i].N()-1;
    int k = int(_rnd.Rannyu() * d);
    vector <gene> c1 = _pop[i].Cut(uno,k,d);
    vector <gene> c2 = _pop[j].Cut(due,k,d);
    
    Sort(c1,due);
    Sort(c2,uno);
    
    _pop[l].Rebuild(uno,c1,k);
    _pop[m].Rebuild(due,c2,k);
    
    _pop[i].Rebuild(uno);
    _pop[j].Rebuild(due);
};

void population::Mutate(){
    for(int i=0;i<_n;i++){
        int m = _rnd.Rannyu() * _pop[i].N();
        double p = _rnd.Rannyu();
        if(p < 0.1) RndSwap(i);
        else if(p < 0.2) RndSwap(i,m/2);
        else if(p < 0.3) RndInv(i);
        else if(p < 0.4) RndShift(i);
    }
};
void population::Evolve(int n){
    Sort();
    for(int k=n;k<_n-1;k+=2){
        Xover(k,k+1);
    }
};
void population::Sort(vector <gene>& guy, vector <gene>& dad){
    std::sort(guy.begin(),guy.end(),[dad](gene& a, gene& b){
        return a.Search(dad) < b.Search(dad);
    });
};
void population::Sort(){
    std::sort(_pop.begin(),_pop.end(),[](const element& a, const element& b){
        return a.L2() < b.L2();
    });
};
void population::Print(int j,int i){
    ofstream fout(_out+"PATH/element_"+to_string(j)+"_"+to_string(i)+".txt");
    for(int j=0;j<_pop[i].N();j++){
        fout<<_pop[i](j).Id();
        for(int k=0;k<_pop[i](j).Pos().n_elem;k++){
            fout<<"\t"<<_pop[i](j).Pos()(k);
        }
        fout<<endl;
    }
    fout.close();
};
void population::L2(int i){
    ofstream fout(_out+"LOSS/l2_gen"+to_string(i)+".txt");
    for(int j=0;j<_n;j++){
        fout<<j<<"\t"<<_pop[j].L2()<<endl;
    }
    fout.close();
};
