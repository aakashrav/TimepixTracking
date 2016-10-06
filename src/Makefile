SOURCES = main.cpp
OBJECTS = $(SOURCES:.c=.o)
INCLUDES = lattice_dbscan.h
CC = g++
CFLAGS = -std=c++11 -Wall -I/Users/AakashRavi/Desktop/Aakash/Education/ROOT/include
LDFLAGS = -L/Users/AakashRavi/Desktop/Aakash/Education/ROOT/lib -lTree -lRIO -lCore \
-lMathCore -lHist -lGraf -lGpad

make: sample_cluster

sample_cluster: $(SOURCES) $(INCLUDES)
	$(CC) $(SOURCES) $(CFLAGS) $(LDFLAGS) -o main

clean:
	rm main