CFLAGS = -Xcompiler -fopenmp -arch compute_20 
CC = nvcc
OBJS = Mesh.cu AUSM_2D.cu Convflux.cu Difflux.cu Pressure.cu Starter.cu Visual.cu
LIBS = -lgomp

ausmPlus_step : ${OBJS}
	${CC} ${CFLAGS} ${OBJS} -o $@ ${LIBS} 

clean:
	-rm -f *.o *.core ausmPlus_step *.csv
