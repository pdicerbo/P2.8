CC = mpic++ #g++
SOURCE=simplemd.cpp Vector.cpp Random.cpp
CXXFLAGS=

all: simplemd.x

simplemd.x: $(SOURCE)
	$(CC) $(SOURCE) $(CXXFLAGS) -O -o simplemd.x

clean:
	rm -f simplemd.x *.o
