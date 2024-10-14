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
#include <algorithm>
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
    int Search(vector <gene> dad);
    
protected:
    string _name;
    int _label;
    vec _pos;
};

class element{
public:
    element(){;};
    element(int dim){
        gene x;
        _dim = dim;
        for(int i=0;i<dim;i++) {
            _dna.push_back(x);
            _sum += (i + 1);
        }
    };
    element(gene* dna,int dim){
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
    void Swap(int i,int j,int m=1); //swap of two m-genes-long fragments
    void Inv(int i, int j);
    
protected:
    vector <gene> _dna;
    int _dim;
    int _sum;
    double _l2;
};

class population{
public:
    population(int N,Random& ran){
        _rnd = ran;
        _n = N;
    }
    population(element pop,int N,Random& ran){
        _rnd = ran;
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
    void Check();
    void RndSwap(int k,int m=1); // randomic swap of two m-genes-long fragments of the k-th element
    void RndInv(int k);
    void Xover(int i, int j);
    //int Search(gene& x, vector <gene> dad);
    void Sort(vector <gene>& guy, vector <gene>& dad);
    void Sort();
    void Print(int i=0);
    void L2(int i=0);
    
protected:
    Random _rnd;
    vector <element> _pop;
    int _n;
};
#endif /* population.h */


