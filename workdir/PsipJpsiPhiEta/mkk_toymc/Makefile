PROGRAM         := toy_mc.exe
#------------------------------------------------------------------------------
# ROOT
#------------------------------------------------------------------------------
ROOTCFLAGS      := $(shell root-config --cflags)
ROOTLDIR        := $(shell root-config --libdir)
ROOTLIBS        := -L$(ROOTLDIR) -lMathMore -lUnuran
ROOTLIBS        += $(shell root-config --libs)
# ROOTGLIBS       := $(shell root-config --glibs)
#------------------------------------------------------------------------------
# GSL
#------------------------------------------------------------------------------
# GSLCFLAGS       := $(shell gsl-config --cflags)
GSLLIBS         := $(shell gsl-config --libs)
#------------------------------------------------------------------------------
CXX             := g++
CXXFLAGS        := -O2 -Wall -pipe -DBATCH
#  CXXFLAGS        := -g -Wall -pipe -DBATCH
#------------------------------------------------------------------------------
# INCLUDES        := -I./
# CXXFLAGS        += $(INCLUDES)
CXXFLAGS        += $(ROOTCFLAGS)
#------------------------------------------------------------------------------
LIBS            := $(ROOTLIBS)
LIBS            += $(GSLLIBS)
#------------------------------------------------------------------------------
SRCS          := ToyMC.cc
#  SRCS          := ToyMC_12.cc
OBJS          := $(SRCS:.cc=.o)

#------------------------------------------------------------------------------
%.o : %.cc
		$(CXX) -c $(CXXFLAGS) $<
#------------------------------------------------------------------------------

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS)
		@echo "Linking $(PROGRAM) ..."
		$(CXX) $(CXXFLAGS) \
		$(OBJS) $(LIBS)\
		-o $(PROGRAM)
		@echo "done"

clean:
		@rm -f *.o $(PROGRAM) core *~

distclean:      clean
		@rm -f make.depends

depend:
		@echo "Generating make.depends"
		@rm -f make.depends
		g++ -MM $(CXXFLAGS) $(SRCS) > make.depends
	        @echo "done"

include $(wildcard *.depends)
