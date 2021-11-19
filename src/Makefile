CC = g++ 
SERIAL = serial
OPENMAP = OpenMp
BIN = bin
OBJ = obj
CFLAGS =  -std=c++11 -fvisibility=hidden -lpthread -O3 -fopenmp -fexceptions

all: $(SERIAL) $(OPENMPFLOYD) $(OPENMP2)

$(SERIAL) : $(OBJ)/serial_vertex_cover.o 
	$(CC) $(CFLAGS) $(OBJ)/serial_vertex_cover.o -o $(BIN)/$(SERIAL)

$(OPENMAP) : $(OBJ)/openmp_vertex_cover.o 
	$(CC) $(CFLAGS) $(OBJ)/openmp_serial_vertex_cover.o -o $(BIN)/$(OPENMAP)

$(OBJ)/serial_vertex_cover.o : ./serial_vertex_cover.cpp ./serial_vertex_cover.h
	$(CC) $(CFLAGS) -c ./serial_vertex_cover.cpp -o $(OBJ)/serial_vertex_cover.o

$(OBJ)/openmp_serial_vertex_cover.o : ./openmp_serial_vertex_cover.cpp ./openmp_serial_vertex_cover.h
	$(CC) $(CFLAGS) -c ./serial_vertex_cover.cpp -o $(OBJ)/serial_vertex_cover.o

clean:
	rm -f $(BIN)/*
	rm -rf $(OBJ)/*.o