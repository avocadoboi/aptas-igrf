COMPILER_FLAGS := -std=c99 -Werror -Wall -Wextra -Iinclude 

build/test: build/aptas-igrf.o test/test.c | build
	gcc $(COMPILER_FLAGS) -o build/test test/test.c build/aptas-igrf.o
build/aptas-igrf.o: include/aptas-igrf.h src/aptas-igrf.c | build
	gcc $(COMPILER_FLAGS) -c -o build/aptas-igrf.o src/aptas-igrf.c
build:
	mkdir build

