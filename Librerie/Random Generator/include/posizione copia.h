//
//  posizione.h
//  
//
//  Created by Giovanni Galafassi on 27/10/21.
//
#ifndef __posizione_h_
#define __posizione_h_
#include <iostream>
#include <cmath>

#include <stdio.h>

class posizione{
public:
    posizione();
    posizione(double, double, double);
    ~posizione();
    
    double GetX()const;
    double GetY()const;
    double GetZ()const;
    double GetR()const;
    double GetPh()const;
    double GetTh()const;
    double GetRh()const;
    void Polar(double rho, double theta, double phi);
    void Set(int i, double val);
    double distance(const posizione&)const;
    posizione operator+(const posizione& v);
    posizione& operator=(const posizione& other);
    
protected:
    double m_x,m_y,m_z;
};

#endif /* posizione_hpp */
