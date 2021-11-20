CC = g++ 
SIMPLE = vanila
SERIAL = serial
OPENMAP = openmp
BIN = bin
OBJ = obj
CFLAGS = -std=c++11 -o3 -fopenmp

all: $(SIMPLE) $(SERIAL) $(OPENMAP)

$(SIMPLE) : $(OBJ)/vanila_vertex_cover.o 
	$(CC) $(CFLAGS) $(OBJ)/vanila_vertex_cover.o -o $(BIN)/$(SIMPLE)

$(SERIAL) : $(OBJ)/serial_vertex_cover.o 
	$(CC) $(CFLAGS) $(OBJ)/serial_vertex_cover.o -o $(BIN)/$(SERIAL)

$(OPENMAP) : $(OBJ)/parallel_vertex_cover.o 
	$(CC) $(CFLAGS) $(OBJ)/parallel_vertex_cover.o -o $(BIN)/$(OPENMAP)

$(OBJ)/vanila_vertex_cover.o : ./vanila_vertex_cover.cpp ./vanila_vertex_cover.h
	$(CC) $(CFLAGS) -c ./vanila_vertex_cover.cpp -o $(OBJ)/vanila_vertex_cover.o

$(OBJ)/serial_vertex_cover.o : ./serial_vertex_cover.cpp ./serial_vertex_cover.h
	$(CC) $(CFLAGS) -c ./serial_vertex_cover.cpp -o $(OBJ)/serial_vertex_cover.o

$(OBJ)/parallel_vertex_cover.o : ./parallel_vertex_cover.cpp ./parallel_vertex_cover.h
	$(CC) $(CFLAGS) -c ./parallel_vertex_cover.cpp -o $(OBJ)/parallel_vertex_cover.o

clean:
	rm -f $(BIN)/*
	rm -rf $(OBJ)/*.o
