//
//  population.cpp
//
//
//  Created by Giovanni Galafassi on 05/10/24.
//
#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdio.h>
#include <map>

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
int gene::Search(vec dad){
    int index = -1;
    for(int i=0;i<dad.size();i++){
        if(_label == int(dad(i))) index = i;
    }
    return index;
};



int& element::GetG(int i){
    return _dna[i];
};
int element::N(){
    return _dim;
};
bool element::Check(){
    int count = 0;
    for(int i=0;i<_dim;i++){
        for(int j=i+1;j<_dim;j++){
            if(_dna[i] == _dna[j]) return false;
        }
        count += _dna[i];
    }
    if (count == _sum) return true;
    else return false;
};
void element::Distances(){
    _d.set_size(_dim,_dim);
    for(int i=0;i<_dim;i++){
        for(int j=0;j<_dim;j++) _d(i,j) = norm(_mp[i] - _mp[j])
    }
};
double element::L2(){
    double l2 = 0.;
    for(int i=0;i<_dim-2;i++){
        l2 += _d(_dna[i],_dna[i+1]);
    }
    l2 += _d(_dna[_dim -1], _dna[0]);
    return l2;
};
int element::PBC(int i,int dim){
    return (i % dim + dim) % dim;
};
vector <int> element::Cut(int i, int j){
    if(j == -1) j = _dim;
    vector <int> app(j - i);
    for(int k=i;k<j;k++) app[k-i] = _dna[k];
    return app;
};
vector <int> element::Cut(vector <int> a,int i, int j){
    if(j == -1) j = _dim;
    vector <int> app(j - i);
    for(int k=i;k<j;k++) app[k-i] = a[PBC(k,a.size())];
    return app;
};
void element::Rebuild(vector <int> r,int i){
    for(int k=i;k<_dim;k++) _dna[k] = r[k-i];
};
void element::Rebuild(vector <int> zip,vector <int> r,int i){
    for(int k=i;k<zip.size();k++) zip[k] = r[k-i];
};
void element::Rebuild(vec labs){
    if(labs.size() != _dim) cerr<<"Wrong indexes' vector size"<<endl;
    else{
        sort(_dna.begin(),_dna.end(),[labs](gene&a, gene& b){
            return a.Search(labs) < b.Search(labs);
        });
    }
    Check();
};
void element::Shift(int n, int m, int i){
    i--;
    vector <int> rap = Cut(1);
    vector <int> uno = Cut(rap,i,i+m);
    vector <int> due = Cut(rap,i+m,i+m+n);
    
    for(int k=i;k<i+n;k++) rap[PBC(k,_dim-1)] = due[k-i];
    for(int k=i+n;k<i+n+m;k++) rap[PBC(k,_dim-1)] = uno[k-i-n];
    Rebuild(rap);
};
void element::Swap(int i,int j,int m){
    i--; j--;
    vector <int> rap = Cut(1);
    vector <int> uno = Cut(rap,i,i+m);
    vector <int> due = Cut(rap,j,j+m);
    
    for(int k=i;k<i+m;k++) rap[PBC(k,_dim-1)] = due[k-i];
    for(int k=j;k<j+m;k++) rap[PBC(k,_dim-1)] = uno[k-j];
    Rebuild(rap);
};
void element::Inv(int i, int j){
    i--;j--;
    vector <int> rap = Cut(1);
    vector <int> uno = Cut(rap,i,j);
    
    for(int k=i;k<j;k++) rap[PBC(k,_dim-1)] = uno[j-k-1];
    Rebuild(rap);
};



element& population::GetEl(int i){
    return _pop[i];
};
element& population::RanEl(){
    return _pop[int(pow(_rnd.Rannyu(),7) * _n)];
};
Random& population::Rnd(){
    return _rnd;
};
void population::Circle(int length){
    element adamo(length);
    for(int i=0;i<length;i++){
        double b = _rnd.Rannyu() * 2. * M_PI;
        vec pos = {sin(b),cos(b)};
        adamo(i) = i;
        adamo[i] = pos;
    }
    adamo.Distances();
    for(int j=0;j<_n;j++) _pop.push_back(adamo);
    Spread();
};
void population::Square(int length){
    element adamo(length);
    for(int i=0;i<length;i++){
        double x = _rnd.Rannyu();
        double y = _rnd.Rannyu();
        vec pos = {x,y};
        adamo(i) = i;
        adamo[i] = pos;
    }
    adamo.Distances();
    for(int j=0;j<_n;j++) _pop.push_back(adamo);
    Spread();
};
void population::File(string file){
    ifstream fin(file);
    string line;
    int length = 0;
    while(getline(fin,line)){
        length++;
    }
    fin.close();
    fin.open(file);
    element adamo(length);
    double x,y;
    for(int i=0;i<length;i++){
        fin >> x >> y;
        vec pos = {x,y};
        adamo(i) = i;
        adamo[i] = pos;
    }
    adamo.Distances();
    for(int j=0;j<_n;j++) _pop.push_back(adamo);
    Spread();
    
};
void population::Spread(){
    for(int j=1;j<_n;j++){
        for(int k=0;k<3;k++){
            RndSwap(j);
            RndSwap(j,2);
            RndInv(j);
            RndShift(j);
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
    vector <int> uno = _pop[i].Cut(1);
    vector <int> due = _pop[j].Cut(1);
    int d = _pop[i].N()-1;
    int k = int(_rnd.Rannyu() * d);
    vector <int> c1 = _pop[i].Cut(uno,k,d);
    vector <int> c2 = _pop[j].Cut(due,k,d);
    
    Sort(c1,due);
    Sort(c2,uno);
    
    _pop[i].Rebuild(c1,k+1);
    _pop[j].Rebuild(c2,k+1);
};
void population::Select(){
    Sort();
    vector <element> app;
    for(int i=0;i<_n;i++)app.push_back(RanEl());
    _pop = app;
};
void population::Mutate(){
    for(int i=0;i<_n;i++){
        double p = _rnd.Rannyu();
        if(p < p_sw) RndSwap(i);
        p = _rnd.Rannyu();
        if(p < p_swm) RndSwap(i,_rnd.Rannyu() * _pop[i].N()/2);
        p = _rnd.Rannyu();
        if(p < p_inv) RndInv(i);
        p = _rnd.Rannyu();
        if(p < p_shf) RndShift(i);
    }
};
void population::Evolve(){
    for(int k=0;k<_n;k++){
        if(_rnd.Rannyu() < p_x) Xover(k,int(_rnd.Rannyu()*_n));
    }
};
void population::Migration(int rank, int cores, int n) {
    
    int dim = _pop[0].N();
    mat sends(dim,n);
    mat recvs(dim,cores*n);

    for(int i=0;i<dim;i++) {
        for(int j=0;j<n;j++) sends(i,j) = _pop[j](i).Id();
    }
    MPI_Gather(sends.memptr(), dim * n, MPI_DOUBLE, recvs.memptr(), dim * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if(rank == 0) recvs = shuffle(recvs,1);

    MPI_Bcast(recvs.memptr(), dim * cores * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for(int j=0;j<n;j++){
        vec a = recvs.col(rank + j);
        _pop[j].Rebuild(a);
    }
};
void population::Sort(vector <gene>& guy, vector <gene>& dad){
    std::sort(guy.begin(),guy.end(),[dad](gene& a, gene& b){
        return a.Search(dad) < b.Search(dad);
    });
};
void population::Sort(){
    std::sort(_pop.begin(),_pop.end(),[](element& a,element& b){
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
