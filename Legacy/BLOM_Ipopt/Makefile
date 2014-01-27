# Copyright (C) 2003, 2010 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# $Id: Makefile.in 2016 2011-06-15 08:53:32Z stefan $

##########################################################################
#    You can modify this example makefile to fit for your own program.   #
#    Usually, you only need to change the five CHANGEME entries below.   #
##########################################################################

# CHANGEME: This should be the name of your executable
EXE = BLOM_NLP

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS =  BLOM_main.o BLOM_NLP_Sparse.o BLOM_NLP.o 

# CHANGEME: Additional libraries
ADDLIBS = Sparse++/lib/libsparse.a Sparse++/lib/libspblas.a \
					Sparse++/lib/libmv.a

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS = -I ../ '-DCOMPLEX=std::complex<double>' -I Sparse++/mv/include 

##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile Ipopt.   #
##########################################################################

# C++ Compiler command
CXX = g++

# C++ Compiler options
CXXFLAGS = -O3 -pipe -DNDEBUG -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wno-unknown-pragmas -Wno-long-long   -DIPOPT_BUILD

# additional C++ Compiler options for linking
CXXLINKFLAGS = 

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
INCL = `PKG_CONFIG_PATH=/Users/JJK/Documents/Research/ipopt/build/lib64/pkgconfig:/Users/JJK/Documents/Research/ipopt/build/lib/pkgconfig:/Users/JJK/Documents/Research/ipopt/build/share/pkgconfig: pkg-config --cflags ipopt` $(ADDINCFLAGS)
#INCL = -I`$(CYGPATH_W) /Users/JJK/Documents/Research/ipopt/build/include/coin`  $(ADDINCFLAGS)

# Linker flags
LIBS = `PKG_CONFIG_PATH=/Users/JJK/Documents/Research/ipopt/build/lib64/pkgconfig:/Users/JJK/Documents/Research/ipopt/build/lib/pkgconfig:/Users/JJK/Documents/Research/ipopt/build/share/pkgconfig: pkg-config --libs ipopt`
##LIBS = -link -libpath:`$(CYGPATH_W) /Users/JJK/Documents/Research/ipopt/build/lib` libipopt.lib -framework vecLib -framework vecLib -lm  -ldl
#LIBS = -L/Users/JJK/Documents/Research/ipopt/build/lib -lipopt `echo -framework vecLib -framework vecLib -lm  -ldl | sed -e 's/ -lmsvcrt / /g'`

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

all: $(EXE)

.SUFFIXES: .cpp .c .o .obj

$(EXE): $(OBJS) $(ADDLIBS)
	bla=;\
	for file in $(OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla $(LIBS) $(ADDLIBS)

$(ADDLIBS):
	cd Sparse++; make sp;

clean:
	rm -rf $(EXE) $(OBJS)
	cd Sparse++; make wipe;

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ $<


.cpp.obj:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `$(CYGPATH_W) '$<'`
