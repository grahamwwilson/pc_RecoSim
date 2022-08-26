# Makefile for the new ntuple processing
# On the CRC cluster, "module load boost" 
main: histset.o recosim.o ParTreeProcessing.C Hungarian.cpp
	g++ -o compiledThreads ParTreeProcessing.C -pthread -I /usr/include/boost `root-config --cflags --libs`

histset.o: histset.C recosim.o Hungarian.h PCTools.h
	g++ -c -pthread  histset.C -I /usr/include/boost `root-config --cflags --libs`

recosim.o: recosim.C recosim.h
	g++ -c -pthread recosim.C `root-config --cflags --libs`

hung.o: Hungarian.cpp Hungarian.h
	g++ -c -pthread Hungarian.cpp `root-config --cflags --libs`
	
test: PairCrossSections.h testCrossSection.cpp
	g++ -o testCrossSection testCrossSection.cpp

clean:
	rm *.o
