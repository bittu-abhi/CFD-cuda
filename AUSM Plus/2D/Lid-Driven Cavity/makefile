CFLAGS = -Xcompiler -fopenmp -arch compute_20
CC = nvcc
OBJS =  compact.cu starter.cu lid.cu sfe.cu visual.cu vte.cu secondW.cu

lidCavity : ${OBJS}
	${CC} ${CFLAGS} ${OBJS} -o $@

clean:
	-rm -f *.o *.core lidCavity *.csv
