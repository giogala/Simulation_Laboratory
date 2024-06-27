//
//  coseno.h
//
//
//  Created by Giovanni Galafassi on 18/11/21.
//
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include "funzione.h"
//#include "seno.h"
#ifndef __coseno_h_
#define __coseno_h_

class coseno: public funzione{
public:
    coseno();
    coseno(double,double,double);
    coseno(const coseno&);
    ~coseno(){;};
    
    double eval(double)const;
    void seta(double);
    void setb(double);
    void setc(double);
    double geta()const;
    double getb()const;
    double getc()const;
    //double xv();
    //double yv();
    seno derivata();
    coseno& operator=(const coseno&);
    coseno& operator=(const seno&);
    
private:
    double m_a, m_b, m_c;
};

#endif /* coseno.h */
