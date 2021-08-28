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
CXXFLAGS        := -g -O2 -Wall -pipe
# CXXFLAGS        := -g -Wall -pipe
#------------------------------------------------------------------------------
# INCLUDES        := -I./
# CXXFLAGS        += $(INCLUDES)
CXXFLAGS        += $(ROOTCFLAGS)
#------------------------------------------------------------------------------
LIBS            := $(ROOTLIBS)
LIBS            += $(GSLLIBS)
#------------------------------------------------------------------------------
SRCS          := toy_mc.cpp
OBJS          := $(SRCS:.cpp=.o)

#------------------------------------------------------------------------------
%.o : %.cpp
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
