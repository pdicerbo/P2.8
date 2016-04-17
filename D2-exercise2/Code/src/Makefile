CC = g++
SOURCE=simplemd.cpp Vector.cpp Random.cpp
CXXFLAGS=

all: simplemd.x

parallel: CXXFLAGS += -D__MPI
parallel: CC = mpic++
parallel: simplemd.x

simplemd.x: $(SOURCE)
	$(CC) $(SOURCE) $(CXXFLAGS) -O3 -o simplemd.x

clean:
	rm -f simplemd.x *.o
