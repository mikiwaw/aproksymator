objects = main.o aproksymator_legendre.o coefficients.o points.o 

legendre: $(objects)
	$(CC) -o legendre $(objects) -lgsl -lgslcblas -lm

main.o : makecof.h
aproksymator_legendre.o : makecof.h
points.o :
coefficients.o :

.PHONY: clean

clean:
	-rm *.o legendre
