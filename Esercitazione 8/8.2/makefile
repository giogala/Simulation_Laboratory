CC = g++
CFLAGS = -Wall -O3 --std=c++11
AFLAGS = -larmadillo


LIB = ../../Librerie/


main : $(LIB)funzione.o $(LIB)random.o $(LIB)datablocking.o $(LIB)metropolis.o $(LIB)library.o main.o
	$(CC) $(CFLAGS) $(AFLAGS) $(LIB)random.o $(LIB)funzione.o $(LIB)datablocking.o $(LIB)metropolis.o $(LIB)library.o main.o -o main

main.o : main.cpp 
	$(CC) $(CFLAGS) -c main.cpp -o main.o

$(LIB)random.o : $(LIB)random.cpp $(LIB)random.h
	$(CC) $(CFLAGS) -c $(LIB)random.cpp -o $(LIB)random.o

$(LIB)funzione.o : $(LIB)datablocking.o $(LIB)funzione.cpp $(LIB)funzione.h
	$(CC) $(CFLAGS) -c $(LIB)funzione.cpp -o $(LIB)funzione.o

$(LIB)datablocking.o : $(LIB)datablocking.cpp $(LIB)datablocking.h
	$(CC) $(CFLAGS) -c $(LIB)datablocking.cpp -o $(LIB)datablocking.o

$(LIB)metropolis.o : $(LIB)metropolis.cpp $(LIB)metropolis.h 
	$(CC) $(CFLAGS) -c $(LIB)metropolis.cpp -o $(LIB)metropolis.o

$(LIB)library.o : $(LIB)library.cpp $(LIB)library.h
	$(CC) $(CFLAGS) -c $(LIB)library.cpp -o $(LIB)library.o

clean :
	rm *.o main $(LIB)*.o

