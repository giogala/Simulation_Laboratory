# Nome dell'eseguibile finale
TARGET = main

# Directory dei sorgenti e delle librerie
SRC_DIR = .
LIB_DIR = ../Librerie

# File sorgenti
SRC_FILES = $(SRC_DIR)/main.cpp
LIB_FILES = $(wildcard $(LIB_DIR)/*.cpp)

# File oggetto
OBJ_FILES = $(SRC_FILES:.cpp=.o) $(LIB_FILES:.cpp=.o)

# Compilatore e flag
CXX = g++
CXXFLAGS = -I$(LIB_DIR) -O3 --std=c++11
AFLAGS = -larmadillo

# Regola predefinita
all: $(TARGET)
	rm $(SRC_DIR)/*.o

# Regola per creare l'eseguibile
$(TARGET): $(OBJ_FILES)
	$(CXX) -o $@ $^ $(AFLAGS)

# Regola per creare i file oggetto
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(CFLAGS)

# Pulizia dei file oggetto e dell'eseguibile
clean:
	rm -f $(OBJ_FILES) $(TARGET)

.PHONY: all clean

