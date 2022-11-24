# Makefile for the new ntuple processing
# # On the CRC cluster, "module load boost" 
main: histsetdata.o datatree.o ParTreeDataProcessing.C Hungarian.cpp
	g++ -o compiledThreadsData ParTreeDataProcessing.C -pthread -I ${BOOSTINCLUDE} `root-config --cflags --libs`

histsetdata.o: histsetdata.C datatree.o Hungarian.h PCToolsData.h PhotonCrossSections.h
	g++ -c -pthread  histsetdata.C -I ${BOOSTINCLUDE} `root-config --cflags --libs`

datatree.o: datatree.C datatree.h
	g++ -c -pthread datatree.C `root-config --cflags --libs`

hung.o: Hungarian.cpp Hungarian.h
	g++ -c -pthread Hungarian.cpp `root-config --cflags --libs`

clean:
	rm *.o
