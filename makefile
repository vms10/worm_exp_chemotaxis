main: main.o WormAgent.o random.o 
	g++ -o main main.o WormAgent.o random.o 
random.o: random.cpp random.h VectorMatrix.h
	g++ -c -O3 random.cpp
WormAgent.o: WormAgent.cpp WormAgent.h random.h VectorMatrix.h
	g++ -c -O3 WormAgent.cpp
main.o: main.cpp WormAgent.h 
	g++ -c -O3 main.cpp
clean: 
	del main.exe *.o 
