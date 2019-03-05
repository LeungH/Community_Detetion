objects=main.o read_network.o maxcliques.o random.o clique_network.o\
	parameter.o community.o caea.o\
	individual.o operator.o evaluation_index.o
LIBS=-lm
CMOEA_SN:$(objects)
	g++ -g -o caeaosn $(objects) $(LIBS)
main.o:main.cpp read_network.h maxcliques.h\
	clique_network.h caea.h
	g++ -g -c main.cpp
read_network.o:read_network.cpp read_network.h
	g++ -g -c read_network.cpp
maxcliques.o:maxcliques.cpp maxcliques.h
	g++ -g -c maxcliques.cpp
random.o:random.cpp random.h
	g++ -g -c random.cpp
clique_network.o:clique_network.cpp clique_network.h clique_network.cpp
	g++ -g -c clique_network.cpp
parameter.o:parameter.cpp 
	g++ -g -c parameter.cpp
community.o:community.cpp community.h
	g++ -g -c community.cpp
individual.o:individual.cpp individual.h
	g++ -g -c individual.cpp
operator.o:operator.cpp individual.h random.h
	g++ -g -c operator.cpp
evaluation_index.o:evaluation_index.cpp evaluation_index.h caea.h
	g++ -g -c evaluation_index.cpp 
caea.o: caea.cpp caea.h individual.h random.h
	g++ -g -c caea.cpp
.PHONY:clean
clean:
	rm caeaosn $(objects)
