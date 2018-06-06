
# location of SuperLU_DIST install
DSuperLUroot = /glade/u/home/klindsay/libs/cheyenne/SuperLU_DIST_5.1.3

# add location of SuperLU_DIST header files to compiler invocation
HEADER = -I$(DSuperLUroot)/SRC

# incorporate SuperLU_DIST make variables, to ensure build is consistent with SuperLU_DIST
include $(DSuperLUroot)/make.inc

# NetCDF flags
NETCDF_CFLAGS = $(shell nc-config --cflags)
NETCDF_LIBS   = $(shell nc-config --libs)

all: bin/gen_A bin/solve_ABglobal bin/solve_ABdist

include Depends

bin/gen_A: gen_A.o $(OBJS)
	$(LOADER) $(LOADOPTS) gen_A.o $(OBJS) $(NETCDF_LIBS) -o $@

bin/solve_ABglobal: solve_ABglobal.o $(OBJS) $(DSUPERLULIB)
	$(LOADER) $(LOADOPTS) solve_ABglobal.o $(OBJS) $(LIBS) $(NETCDF_LIBS) -lm -o $@

bin/solve_ABdist: solve_ABdist.o $(OBJS) $(DSUPERLULIB)
	$(LOADER) $(LOADOPTS) solve_ABdist.o $(OBJS) $(LIBS) $(NETCDF_LIBS) -lm -o $@

.c.o:
	$(CC) $(CFLAGS) $(NETCDF_CFLAGS) $(CDEFS) $(HEADER) -c $< $(VERBOSE)

clean:
	rm -f *.o bin/gen_A bin/solve_ABglobal bin/solve_ABdist
