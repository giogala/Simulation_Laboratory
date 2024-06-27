#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "../../Librerie/Random Generator/random.h"

using namespace std;

int main (int argc, char** argv) {
    if(argc<3){
        cerr<<"Errore: uso del programma "<<argv[0]<<" <Numero_di_lanci> <Numero_di_blocchi"<<endl;
        return -1;
    }
    
    Random rnd;                                 //Imposto il generatore di numeri casuali
    int seed[4]={5,2,3,4};
    rnd.SetRandom(seed,5,10);
    
    int N=atoi(argv[1]);
    int M=atoi(argv[2]);
    double pi=3.;
    double pi2=0.;
    double ti=0.;
    double ti2=0.;
    double L=1.;
    double d=1.5;

    
    //File di output
    ofstream fout("Data.txt");
    fout<<"π"<<"\t"<<"sigma-π"<<endl;
    /*                                                                              p2*
     -------------------------------------------------------------------------------|-----------------
     ^y                                     p2                                      p1*
     |                                     /                     p1——p2
     |                                    p1
     |               p1*                                                                          x
     |-----------------\----------------------------------------------------------------------------->
                        p2*
     Il problema è invariante per traslazioni discrete lungo l'asse x per cui posso ridurmi al problema
     tra solo due righe. Date le due righe distanti d, considero il punto p1, un estremo del bastoncino,
     generato, wlog, nello spazio tra le due righe. Il problema sta nella genazione nella direzione del
     bastoncino. Il problema è invariante per traslazioni lungo l'asse x, ciò significa innanzitutto che
     di p1 ci interessa unicamente l'ordinata. Per quanto riguarda l'angolo, non potendo chiedere che
     venga generato tra 0 e 2π o π, tento un altro approcio: genero angoli casuali tra 0 e la mia stima
     di π e man mano che reitero la correggo. Questo metodo per quanto affascinante non coverge e continua
     ad oscillare attorno al valore atteso e inoltre non mi permette di fornire una stima dell'errore.
     Avrei voluto inoltre fare uso veramente di due soli numeri causuali ma non utilizzare π mi impedisce
     di sfruttare veramente gli angoli. Per determinare la direzione causale genero punti in un quadrato
     e seleziono unicamente quelli con distanza minore di 1, normalizzo, ed ecco un versone distribuito
     isotropicamente.
     
     */
    //Parte 1
    
    for(int j=0;j<M;j++){
        int c=0;                            //contatore di punti del quadrato che stanno dentro al cerchio
        int H=0;                            //contatore degli hit
        for(int i=0;i<N/M;i++){
            //genero un punto in un quadrato di lato 2
            double x = rnd.Rannyu(-1.,1.);
            double y = rnd.Rannyu(-1.,1.);
            double r = x*x+y*y;
            if(r<=1){                                       //se sta dentro il cerchio di raggio 1
                double py= rnd.Rannyu()*d;                  //genero una y casuale in [0,d)
                double y2=py+L*y/sqrt(r);                   //mi sposto della y rinormalizzata
                c++;
                if(y2>d or y2<0) H++;
            }
        }
        pi = (double)(2.*L*c)/(H*d);
        pi2 = pi*pi;
        ti=(ti*j+pi)/(j+1);
        ti2=(ti2*j+pi2)/(j+1);
        fout<<ti<<"\t"<<sqrt((ti2-ti*ti)/j)<<endl;
    }
    
    fout.close();
    cout<<"Buffon π Experiment"<<endl;
    cout<<"------------------------------------------"<<endl;
    cout<<"^y                                     p2 "<<endl;
    cout<<"|                                     /   "<<endl;
    cout<<"|                                    p1   "<<endl;
    cout<<"|                p1*                      "<<endl;
    cout<<"|-----------------|---------------------->"<<endl;
    cout<<"                 p2*                    x "<<endl;
    return 0;
}

