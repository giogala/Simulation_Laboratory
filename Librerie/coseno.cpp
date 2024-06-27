//
//  coseno.cpp
//
//
//  Created by Giovanni Galafassi on 18/11/21.
//

#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
//#include "seno.h"
#include "coseno.h"

using namespace std;

coseno::coseno(){
    m_a=0; m_b=0; m_c=0;
}
coseno::coseno(double a,double b, double c){
    m_a=a; m_b=b; m_c=c;
}
coseno::coseno(const seno& s){
    m_a=s.geta(); m_b=s.getb(); m_c=s.getc();
}


/*parabola::parabola(const parabola& p){
    m_2=p.geta(); m_1=p.getb(); m_0=p.getc();
}*/

double coseno::eval(double x)const {
    return m_a*cos(m_b*x + m_c);
}
void coseno::seta(double v){m_a=v;}
void coseno::setb(double v){m_b=v;}
void coseno::setc(double v){m_c=v;}
double coseno::geta()const{return m_a;}
double coseno::getb()const{return m_b;}
double coseno::getc()const{return m_c;}

//double parabola::xv(){return -m_1/(2*m_2);}
//double parabola::yv(){return eval(xv());}
seno coseno::derivata(){
    coseno c(-geta()*getb(), getb(), getc());
    return c;
}

coseno& coseno::operator=(const coseno& s){
    m_a=s.geta(); m_b=s.getb(); m_c=s.getc();
    return *this;
}
coseno& coseno::operator=(const seno& s){
    m_a=s.geta(); m_b=s.getb(); m_c=s.getc()+M_PI/2.;
    return *this;
}
