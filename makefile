CC = mpicc
CFLAGS= -g -Wall -gdwarf-2 -std=c99 
LDFLAGS= -lm


all:rma_fence rma_pscw

sendrecv:sendrecv.o jacobi.o decomp.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
nonblocking:nonblocking.o jacobi.o decomp.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
rma_fence: rma_fence.o jacobi.o decomp.o 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
rma_pscw: rma_pscw.o jacobi.o decomp.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)





.PHONY:clean

clean:
	rm -f rma_fence.o rma_pscw.o rma_fence rma_pscw jacobi.o decomp.o sendrecv.o sendrecv nonblocking.o nonblocking
