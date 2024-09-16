//
//  posizione.cpp
//  
//
//  Created by Giovanni Galafassi on 27/10/21.
//

#include <iostream>
#include <cmath>
#include "posizione.h"

posizione::posizione(){
    m_x=0;
    m_y=0;
    m_z=0;
}

posizione::posizione(double a, double b, double c){
    m_x=a;
    m_y=b;
    m_z=c;
}

posizione::~posizione(){
}
//coordinate cartesiane
double posizione::GetX()const{
    return m_x;
}
double posizione::GetY()const{
    return m_y;
}
double posizione::GetZ()const{
    return m_z;
}
//coordinate sferiche
double posizione::GetR()const{
    return sqrt(pow(m_x,2)+pow(m_y,2)+pow(m_z,2));
}
double posizione::GetPh()const{
    return atan2(m_y,m_x);
}
double posizione::GetTh()const{
    return acos(m_z/GetR());
}
//coordianta cilindrica
double posizione::GetRh()const{
    return sqrt(m_x*m_x+m_y*m_y);
}
// imposto coordinate sferiche
void posizione::Polar(double rho, double theta, double phi){
    m_x = rho*sin(theta)*cos(phi);
    m_y = rho*sin(theta)*sin(phi);
    m_z = rho*cos(theta);
}
//distanza da un altro punto
double posizione::distance(const posizione& v)const{
    double x=m_x-v.GetX();
    double y=m_y-v.GetY();
    double z=m_z-v.GetZ();
    
    return sqrt(x*x+y*y+z*z);
}
void posizione::Set(int i, double val){
    if(i==0){ m_x=val;}
    if(i==1){ m_y=val;}
    if(i==2){ m_z=val;}
}
// implemento la somma
posizione posizione::operator+(const posizione& v){
    return posizione(m_x+v.GetX(),m_y+v.GetY(),m_z+v.GetZ());
}
// assegnazione
posizione& posizione::operator=(const posizione& other) {
    if (this != &other) {
        m_x = other.GetX();
        m_y = other.GetY();
        m_z = other.GetZ();
    }
    return *this;
}
