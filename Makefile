all: HeuristicDetMot

HeuristicDetMot: src/main.cpp
	g++ src/main.cpp -o HeuristicDetMot

run: HeuristicDetMot
	./HeuristicDetMot