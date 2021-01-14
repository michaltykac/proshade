## /***************************************************************************
##   **************************************************************************
##   
##   SOFT: SO(3) Fourier Transforms
##   Version 2.0
## 
##   Copyright (c) 2003, 2004, 2007 Peter Kostelec, Dan Rockmore
##   
##   This file is part of SOFT.
## 
##   SOFT is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 3 of the License, or
##   (at your option) any later version.
## 
##   SOFT is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
## 
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.
##   
##   See the accompanying LICENSE file for details.
##   
##   ************************************************************************
##   ************************************************************************/
## 

## define WALLCLOCK if want walltime and not cputime
## CFLAGS = -Wall -g -DWALLCLOCK
## CFLAGS = -Wall -g
## CFLAGS = -O3 -pg
## CFLAGS = -O3
## CFLAGS = -Wall -g -I./include $(FFTWINC)
CFLAGS = -Wall -O3 -I./include $(FFTWINC)

FFTWLIB := $(FFTWLIB1) $(FFTWLIB2)

LIBNAME = libsoft1.a
LDFLAGS = -L. -lsoft1 -lm $(FFTWLIB)

COMSRC = $(wildcard common/*.c)
COMOBJ = $(patsubst %.c,%.o,$(COMSRC))

LIBSRC = $(COMSRC) \
	$(wildcard lib1/*.c)
LIBOBJ = $(patsubst %.c,%.o,$(LIBSRC))


EXELIST = test_wigSpec test_genWig test_Wigner_Analysis \
	test_Wigner_Synthesis test_Wigner_Naive test_Wigner_angle
EXELOC = $(patsubst %,bin/%,$(EXELIST))

EXE1LIST = test_soft_fftw test_soft_fftw_for test_soft_fftw_inv \
	test_soft_fftw_pc test_soft_fftw_nt test_s2_rotate_fftw \
	test_s2_rotate_fftw_wrap test_s2_rotate_fftw_mem \
	test_soft_fftw_correlate2 test_soft_fftw_correlate2_wrap
EXE1LOC = $(patsubst %,bin/%,$(EXE1LIST))

###################################################################
##################################################################
######
######              things that can be made
######
##################################################################
##################################################################

# things that can be made

depend:
	@echo $(FFTWINC)
	@echo $(FFTWLIB1)
	@echo $(FFTWLIB2)
	@echo $(FFTWLIB)
	makedepend -I. ${FFTWINC} $(SRC)

clean:
	rm -f {common,lib1,examples,examples1}/*.o

vclean:
	rm -f libsoft1.a $(EXELOC) $(EXE1LOC) {common,lib1,examples,examples1}/*.o

libsoft1.a: $(LIBOBJ)
	ar vru $(LIBNAME) $(LIBOBJ)
	ranlib $(LIBNAME)

lib:	libsoft1.a


#
# making the example Wigner routines
#

examples: $(EXELIST)

#
# making the example SO(3) FT routines
#

examples1: $(EXE1LIST)


tests: examples examples1


all: lib examples examples1


########################################################
#
# making the executables
#
########################################################

#
# discrete Wigner transform routines
#

test_wigSpec: examples/test_wigSpec.o libsoft1.a
	$(CC) $(CFLAGS) examples/test_wigSpec.o \
	$(LDFLAGS) -o bin/test_wigSpec

test_genWig: examples/test_genWig.o libsoft1.a
	$(CC) examples/test_genWig.o \
	$(LDFLAGS) -o bin/test_genWig

test_Wigner_angle: examples/test_Wigner_angle.o libsoft1.a
	$(CC) $(CFLAGS) examples/test_Wigner_angle.o \
	$(LDFLAGS) -o bin/test_Wigner_angle

test_Wigner_Analysis: examples/test_Wigner_Analysis.o libsoft1.a
	$(CC) $(CFLAGS) examples/test_Wigner_Analysis.o \
	$(LDFLAGS) -o bin/test_Wigner_Analysis

test_Wigner_Synthesis: examples/test_Wigner_Synthesis.o libsoft1.a
	$(CC) $(CFLAGS) examples/test_Wigner_Synthesis.o \
	$(LDFLAGS) -o bin/test_Wigner_Synthesis

test_Wigner_Naive: examples/test_Wigner_Naive.o libsoft1.a
	$(CC) $(CFLAGS) examples/test_Wigner_Naive.o \
	$(LDFLAGS) -o bin/test_Wigner_Naive


#
# SO(3) routines that depend on FFTW
#


test_soft_fftw: examples1/test_soft_fftw.o libsoft1.a
	$(CC) $(CFLAGS) examples1/test_soft_fftw.o \
	$(LDFLAGS) ${FFTWLIB} -o bin/test_soft_fftw

test_soft_fftw_for: examples1/test_soft_fftw_for.o libsoft1.a
	$(CC) $(CFLAGS) examples1/test_soft_fftw_for.o \
	$(LDFLAGS) ${FFTWLIB} -o bin/test_soft_fftw_for

test_soft_fftw_inv: examples1/test_soft_fftw_inv.o libsoft1.a
	$(CC) $(CFLAGS) examples1/test_soft_fftw_inv.o \
	$(LDFLAGS) ${FFTWLIB} -o bin/test_soft_fftw_inv

test_soft_fftw_pc: examples1/test_soft_fftw_pc.o libsoft1.a
	$(CC) $(CFLAGS) examples1/test_soft_fftw_pc.o \
	$(LDFLAGS) ${FFTWLIB} -o bin/test_soft_fftw_pc

test_soft_fftw_nt: examples1/test_soft_fftw_nt.o libsoft1.a
	$(CC) $(CFLAGS) examples1/test_soft_fftw_nt.o \
	$(LDFLAGS) ${FFTWLIB} -o bin/test_soft_fftw_nt

test_s2_rotate_fftw: examples1/test_s2_rotate_fftw.o libsoft1.a
	$(CC) $(CFLAGS) examples1/test_s2_rotate_fftw.o \
	$(LDFLAGS) ${FFTWLIB} -o bin/test_s2_rotate_fftw

test_s2_rotate_fftw_wrap: examples1/test_s2_rotate_fftw_wrap.o libsoft1.a
	$(CC) $(CFLAGS) examples1/test_s2_rotate_fftw_wrap.o \
	$(LDFLAGS) ${FFTWLIB} -o bin/test_s2_rotate_fftw_wrap

test_s2_rotate_fftw_mem: examples1/test_s2_rotate_fftw_mem.o libsoft1.a
	$(CC) $(CFLAGS) examples1/test_s2_rotate_fftw_mem.o \
	$(LDFLAGS) ${FFTWLIB} -o bin/test_s2_rotate_fftw_mem

test_soft_fftw_correlate2: examples1/test_soft_fftw_correlate2.o libsoft1.a
	$(CC) $(CFLAGS) examples1/test_soft_fftw_correlate2.o \
	$(LDFLAGS) ${FFTWLIB} -o bin/test_soft_fftw_correlate2

test_soft_fftw_correlate2_wrap: examples1/test_soft_fftw_correlate2_wrap.o libsoft1.a
	$(CC) $(CFLAGS) examples1/test_soft_fftw_correlate2_wrap.o \
	$(LDFLAGS) ${FFTWLIB} -o bin/test_soft_fftw_correlate2_wrap



# and now for LOTS OF dependencies ...

# DO NOT DELETE THIS LINE -- make depend depends on it.
