all:
	g++ main.cpp -o sph_tree.bin -O3 -g -lglut -lGLU -funroll-loops -ffast-math -malign-double -lpng
