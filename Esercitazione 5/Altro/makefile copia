#LIB_DIR := ../Librerie/Random\ Generator/
#INC_DIR := ../Librerie/Random\ Generator/include/

# Utilizzo di queste variabili per definire le variabili LIBS e INCS
#LIBS := -L$(LIB_DIR)
#INCS := -I$(INC_DIR)
CFLAGS = -O3 --std=c++11
AFLAGS = -larmadillo

rmo : main
	rm *.o

main : main.o random.o funzione.o metropolis.o
	g++ -o main main.o random.o $(AFLAGS)

main.o : main.cpp
	g++ -c main.cpp -o main.o $(CFLAGS)

random.o : random.cpp random.h
	g++ -c random.cpp -o random.o $(CFLAGS)
	
funzione.o : funzione.h
	g++ -c funzione.h -o funzione.o $(CFLAGS)

metropolis.o : metropolis.h
	g++ -c metropolis.h -o metropolis.o $(CFLAGS)

#$(INC_DIR)random.h : ../Librerie/Random\ Generator/random.h
#	mkdir -p $(INC_DIR)
#	cp ../Librerie/Random\ Generator/random.h $(INC_DIR)
	
#$(INC_DIR)posizione.h : ../Librerie/Random\ Generator/posizione.h
#	mkdir -p $(INC_DIR)
#	cp ../Librerie/Random\ Generator/posizione.h $(INC_DIR)

clean:
	rm -f main



