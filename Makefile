COMPILER_FLAGS := -std=c99 -Wall -Wextra -Iinclude -Iold-implementation -O2
CC := gcc
# COMPILER_FLAGS := -DCPP_ASSOC_LEGENDRE_TESTING -std=c++23 -Wall -Wextra -Iinclude 
# CC := g++

build/test: build/aptas-igrf.o build/geomag70.o test/test.c | build
	$(CC) $(COMPILER_FLAGS) -o $@ test/test.c build/aptas-igrf.o build/geomag70.o -lm -O2
build/aptas-igrf.o: include/aptas-igrf.h src/aptas-igrf.c | build
	$(CC) $(COMPILER_FLAGS) -c -o $@ src/aptas-igrf.c -O2
build/geomag70.o: old-implementation/geomag70.h old-implementation/geomag70.c | build 
	$(CC) $(COMPILER_FLAGS) -c -o $@ old-implementation/geomag70.c -O2
	
build:
	mkdir build
clean:
	rm -r build
