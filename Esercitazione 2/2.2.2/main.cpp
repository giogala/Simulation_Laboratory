
#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "../../Librerie/Random Generator/random.h"
#include "../../Librerie/Random Generator/posizione.h"

using namespace std;

int main (int argc, char** argv) {
    if(argc<3){
        cerr<<"Errore: uso del programma "<<argv[0]<<" <Numero_di_passi> <Numero_di_repliche> <Numero_blocchi_di_repliche"<<endl;
        return -1;
    }
    
    Random rnd;                                 //Imposto il generatore di numeri casuali
    int seed[4]={5,2,3,4};
    rnd.SetRandom(seed,5,100);
    
    //File di output
    ofstream fout("Data.txt");
    ofstream fout2("Data2.txt");
    fout<<"$<R^2>_{RW}$ "<<"\t"<<"devstd"<<endl;
    fout2<<"X"<<"\t"<<"Y"<<"\t"<<"Z"<<endl;
    
    int N=atoi(argv[1]);                        //Numero di passi
    int M=atoi(argv[2]);                        //Numero di repliche
    int L=atoi(argv[3]);                        //Numero di blocchi in cui suddividerwe le M repliche
    int O=M/L;                                  //Numero di repliche per blocco
    
    // Vettore con tutte le posizioni di tutte le repliche del RW inizializzate nell'origine
    //vector <vector <double> > pos;
    vector <posizione> pos;
    
    double r2=0;            //variabile locale nel blocco per il calcolo del valor medio
    double r4=0;            //variabile locale nel blocco per il calcolo della varianza
    double t2=0;            //variabile globale per il calcolo del valor medio dei valori medi a i fissato
    double t4=0;            //variabile globale per il calcolo della dev. std della media a i fissato
    
    //vector<double> r(3,0.);
    posizione r;
    for(int i=0;i<M;i++){
        pos.push_back(r);
    }
    
    // iterazione sul numero di passi
    for(int i=0;i<N;i++){
        // iterazione sui blocchi in cui ho diviso le repliche
        for(int k=0;k<L;k++){
            //iterazione sulle repliche all'interno di un blocco
            for(int j=O*k;j<O*(k+1);j++){
                double th=rnd.Rannyu(0,M_PI);           //genero casualmente theta in [0,π)
                double phi=rnd.Rannyu(0,2*M_PI);        //genero casualmente phi in [0,2π)
                //posizione r();
                r.Polar(1.,th,phi);
                pos[j]=pos[j]+r;                        //evolvo posizione della replica
                
                
                // stampo su Data2 l'evoluzione della prima replica per visualizzarla
                if(j==0){
                    fout2<<pos[j].GetX()<<"\t"<<pos[j].GetY()<<"\t"<<pos[j].GetZ()<<endl;
                }
                //incremento il contatore locale di r2 con la norma del vettore posizione della j-esima replica fratto il numero di repliche nel blocco
                
                r2+=pow(pos[j].GetR(),2)/O;
                
            }
            r2=sqrt(r2);
            r4=r2*r2;
            t2+=r2/L;
            t4+=r4/L;
            r2=0; r4=0;
        }
        double e = sqrt((t4-t2*t2)/(L-1));      // il "/L" è dovuto al fatto che sto mediando su L blocchi
        fout<<t2<<"\t"<<e<<endl;
        t2=0; t4=0;
    }
    
    fout.close();
    fout2.close();
    return 0;
}

