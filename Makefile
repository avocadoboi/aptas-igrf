COMPILER_FLAGS := -std=c99 -Wall -Wextra -Iinclude
# COMPILER_FLAGS := -std=c++23 -Wall -Wextra -Iinclude 
CC := gcc

build/test: build/aptas-igrf.o test/test.c | build
	$(CC) $(COMPILER_FLAGS) -o $@ test/test.c build/aptas-igrf.o -lm
build/aptas-igrf.o: include/aptas-igrf.h src/aptas-igrf.c | build
	$(CC) $(COMPILER_FLAGS) -c -o $@ src/aptas-igrf.c
build:
	mkdir build
clean:
	rm -r build
