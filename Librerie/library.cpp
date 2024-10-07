#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>  // Per std::fflush
#include "library.h"

using namespace std;

void Progress_Bar(int progress, int total, int bar_width) {
    // Calcola la percentuale di completamento
    float percentage = static_cast<float>(progress) / total;
    int pos = static_cast<int>(bar_width * percentage);

    // Crea la barra di progresso
    string bar;
    for (int i = 0; i < bar_width; ++i) {
        if (i <= pos)
            bar += "◼︎";  // Puoi sostituire con '=' o altri caratteri compatibili
        else
            bar += " ";
    }

    // Stampa la barra di progresso e la percentuale
    cout << "|" << bar << "| " << int(percentage * 100) << " %\r";
    cout.flush();  // Assicurati che la barra venga stampata immediatamente
}

double SetProp(string file,string prop) {
    ifstream fin(file);
    string input;
    double r = 0.;
    
    if(fin.is_open()){
        while(!fin.eof()){
            fin >> input;
            if( input == prop ) fin >> r;
        }
    } else cerr<<"Can't open "<<file<<endl;
    fin.close();
    return r;
}
