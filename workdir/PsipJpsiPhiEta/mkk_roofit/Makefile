ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --glibs)
ROOTLIBS     += -lRooFitCore -lRooFit -lMinuit
#---------------------------------------------------
CXX          := g++
# CXXFLAGS     := -O2 -Wall -pipe
CXXFLAGS     = -g -Wall -pipe
#---------------------------------------------------
LD           := gcc
#------------------------------------------------------------------------------
INCLUDES     := -I./
CXXFLAGS     += $(INCLUDES)
CXXFLAGS     += $(ROOTCFLAGS)
#------------------------------------------------------------------------------
LIBS         := $(ROOTLIBS)
#------------------------------------------------------------------------------
SRCS         := $(wildcard *.cxx)
OBJS         := $(SRCS:.cxx=.o)
PROGRAM      := mkk_fit
#------------------------------------------------------------------------------
# rules
%.o : %.cxx
		$(CXX) -c $(CXXFLAGS) $<
#------------------------------------------------------------------------------

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS)
		@echo "Linking $(PROGRAM) ..."
		$(CXX) $(CXXFLAGS) $(OBJS) \
		$(LIBS) \
		-o $(PROGRAM)
		@echo "done"

clean:
		@rm -f *.o

distclean:      clean
		@rm -f $(PROGRAM) 

