LIB_DIR := ../../Librerie/Random\ Generator/
INC_DIR := ../../Librerie/Random\ Generator/include/

# Utilizzo di queste variabili per definire le variabili LIBS e INCS
LIBS := -L$(LIB_DIR)
INCS := -I$(INC_DIR)

main : main.o random.o
	g++ -o main main.o random.o $(LIBS)

main.o : main.cpp
	g++ -c main.cpp -o main.o $(INCS)

random.o : $(LIB_DIR)random.cpp $(INC_DIR)random.h
	g++ -c $(LIB_DIR)random.cpp -o random.o $(INCS)

$(INC_DIR)random.h : ../../Librerie/Random\ Generator/random.h
	mkdir -p $(INC_DIR)
	cp ../../Librerie/Random\ Generator/random.h $(INC_DIR)

clean:
	rm -f main main.o random.o



