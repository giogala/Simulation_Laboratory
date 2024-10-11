#ifndef __library_h_
#define __library_h_
#include <iostream>
#include <string>
#include <cstdio> 
#include <fstream>
using namespace std;


void Progress_Bar(int progress, int total, int bar_width = 50);
double SetProp(string file,string prop);
int fact(int i);

#endif /* library.h */
