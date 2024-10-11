//
//  population.h
//
//
//  Created by Giovanni Galafassi on 09/05/24.
//
#define _USE_MATH_DEFINES
#ifndef __population_h_
#define __population_h_
#include <iostream>
#include <cmath>
#include <armadillo>
#include <stdio.h>
#include <vector>
#include "library.h"
#include "random.h"

using namespace std;
using namespace arma;

class gene{
public:
    gene(){;};
    gene(string name,int i,vec x){
        _name = name;
        _label = i;
        _pos = x;
    };
    ~gene(){;};
    
    void Naming(string n);
    void SetId(int i);
    void Placing(vec x);
    void Init(string name,int i,vec x);
    int Id();
    vec Pos()const;
    
protected:
    string _name;
    int _label;
    vec _pos;
};

class element{
public:
    element(){;};
    element(int dim,Random& ran){
        _rnd = ran;
        gene x;
        _dim = dim;
        for(int i=0;i<dim;i++) {
            _dna.push_back(x);
            _sum += (i + 1);
        }
    };
    element(gene* dna,int dim,string rnd_dir){
        _rnd.initRnd(rnd_dir);
        _dim = dim;
        for(int i=0;i<dim;i++) {
            _dna.push_back(dna[i]);
            _sum += (i + 1);
        }
    };
    ~element(){;};
    
    gene& operator()(int i) {return GetG(i);};
    gene& GetG(int i);
    int N();
    bool Check();
    double L2()const;
    void SetL2();
    int PBC(int i,int dim);
    
    vector <gene> Cut(int i, int j=-1);
    vector <gene> Cut(vector <gene> a,int i, int j=-1);
    void Rebuild(vector <gene> r,int i = 1);
    void Shift(int n, int m, int i); // +n shift of m contiguos genes from i-th position on
    void RndSwap(); // randomic swap of two genes
    void RndSwap(int m); // randomic swap of two m-genes-long fragments
    void Inv(int i, int j);
    void RndInv();
    
protected:
    Random _rnd;
    vector <gene> _dna;
    int _dim;
    int _sum;
    double _l2;
};

class population{
public:
    population(int N,string rnd_dir){
        _rnd.initRnd(rnd_dir);
        _n = N;
    }
    population(element pop,int N,string rnd_dir){
        _rnd.initRnd(rnd_dir);
        _n = N;
        for(int i=0;i<N;i++) _pop.push_back(pop);
    };
    ~population(){;};
    
    element& operator()(int i) {return GetEl(i);};
    element& operator()() {return RanEl();};
    element& GetEl(int i);
    element& RanEl();
    void Circle(int length);
    void Spread();
    int Check();
    void Xover(int i, int j);
    int Search(gene x, vector <gene> dad);
    void Sort(vector <gene> guy, vector <gene> dad);
    void Sort();
    void Print(int i=0);
    
protected:
    Random _rnd;
    vector <element> _pop;
    int _n;
};
#endif /* population.h */


