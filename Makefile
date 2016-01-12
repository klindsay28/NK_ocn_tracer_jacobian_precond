include make.inc

REGENW := FALSE

HEADER = -I$(DSuperLUroot)/SRC

all: bin/gen_A bin/solve_ABglobal bin/solve_ABdist

include Depends

bin/gen_A: gen_A.o $(OBJS)
	$(LOADER) $(LOADOPTS) gen_A.o $(OBJS) -o $@

bin/solve_ABglobal: solve_ABglobal.o $(OBJS) $(DSUPERLULIB)
	$(LOADER) $(LOADOPTS) solve_ABglobal.o $(OBJS) $(LIBS) -lm -o $@

bin/solve_ABdist: solve_ABdist.o $(OBJS) $(DSUPERLULIB)
	$(LOADER) $(LOADOPTS) solve_ABdist.o $(OBJS) $(LIBS) -lm -o $@

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) $(HEADER) -c $< $(VERBOSE)

clean:
	rm -f *.o bin/gen_A bin/solve_ABglobal bin/solve_ABdist
