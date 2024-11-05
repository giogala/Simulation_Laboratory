# Nome dell'eseguibile finale
TARGET = main

# Directory dei sorgenti e delle librerie
SRC_DIR = .
LIB_DIR = ../Librerie

# File sorgente principale
SRC_FILE = $(SRC_DIR)/main.cpp

# File sorgenti delle librerie specifiche
LIB_SOURCES = $(LIB_DIR)/random.cpp $(LIB_DIR)/library.cpp $(LIB_DIR)/datablocking.cpp $(LIB_DIR)/metropolis.cpp $(LIB_DIR)/funzione.cpp

# File oggetto
OBJ_FILES = $(SRC_FILE:.cpp=.o) $(LIB_SOURCES:.cpp=.o)

# Compilatore e flag
CXX = g++
CXXFLAGS = -I$(LIB_DIR) -O3 --std=c++11
AFLAGS = -larmadillo

# Regola predefinita
all: $(TARGET)

# Regola per creare l'eseguibile
$(TARGET): $(OBJ_FILES)
	$(CXX) -o $@ $^ $(AFLAGS)

# Regola per creare i file oggetto dai file sorgenti
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(LIB_DIR)/%.o: $(LIB_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Pulizia dei file oggetto e dell'eseguibile
clean:
	rm -f $(OBJ_FILES) $(TARGET)

.PHONY: all clean
