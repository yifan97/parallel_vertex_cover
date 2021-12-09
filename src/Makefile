CC = g++ 
SIMPLE = vanila
SERIAL = serial
OPENMAP = openmp
MPI = mpi
BIN = bin
OBJ = obj
CFLAGS = -std=c++11 -o3 -fopenmp
CXX = mpic++
CXXFLAGS = -I. -O3

all: $(SIMPLE) $(SERIAL) $(OPENMAP) $(MPI)

$(SIMPLE) : $(OBJ)/vanila_vertex_cover.o 
	$(CC) $(CFLAGS) $(OBJ)/vanila_vertex_cover.o -o $(BIN)/$(SIMPLE)

$(SERIAL) : $(OBJ)/serial_vertex_cover.o 
	$(CC) $(CFLAGS) $(OBJ)/serial_vertex_cover.o -o $(BIN)/$(SERIAL)

$(OPENMAP) : $(OBJ)/parallel_vertex_cover.o 
	$(CC) $(CFLAGS) $(OBJ)/parallel_vertex_cover.o -o $(BIN)/$(OPENMAP)

$(MPI) : $(OBJ)/MPI_vertex_cover.o 
	$(CXX) $(CXXFLAGS) $(OBJ)/MPI_vertex_cover.o -o $(BIN)/$(MPI) 

$(OBJ)/vanila_vertex_cover.o : ./vanila_vertex_cover.cpp ./vanila_vertex_cover.h
	$(CC) $(CFLAGS) -c ./vanila_vertex_cover.cpp -o $(OBJ)/vanila_vertex_cover.o

$(OBJ)/serial_vertex_cover.o : ./serial_vertex_cover.cpp ./serial_vertex_cover.h
	$(CC) $(CFLAGS) -c ./serial_vertex_cover.cpp -o $(OBJ)/serial_vertex_cover.o

$(OBJ)/parallel_vertex_cover.o : ./parallel_vertex_cover.cpp ./parallel_vertex_cover.h
	$(CC) $(CFLAGS) -c ./parallel_vertex_cover.cpp -o $(OBJ)/parallel_vertex_cover.o

$(OBJ)/MPI_vertex_cover.o : ./MPI_vertex_cover.cpp ./MPI_vertex_cover.h 
	$(CXX) $(CXXFLAGS) -c ./MPI_vertex_cover.cpp -o $(OBJ)/MPI_vertex_cover.o 

clean:
	rm -f $(BIN)/*
	rm -rf $(OBJ)/*.o
