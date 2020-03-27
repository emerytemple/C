OBJECTS = main.o geometry.o grid.o cfd.o mymath.o
CFLAGS = -lm -Wall

main: $(OBJECTS)
	cc -o main $(OBJECTS) $(CFLAGS)

geometry.o:
	cc -c geometry.c
grid.o:
	cc -c grid.c
cfd.o:
	cc -c cfd.c
mymath.o:
	cc -c mymath.c

.PHONY: clean
clean:
	-rm main $(OBJECTS)
