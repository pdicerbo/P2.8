CC=mpic++
SOURCE=simplemd.cpp Vector.cpp Random.cpp

all: simplemd.x

simplemd.x: $(SOURCE)
	$(CC) $(SOURCE) -O -o simplemd.x

clean:
	rm -f simplemd.x *.o
