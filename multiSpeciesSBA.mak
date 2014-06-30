.PHONY : clean

#vpath_src=.. ../../randlib/src 
vpath %.c    $(vpath_src)
vpath %.cpp  $(vpath_src)
vpath %.hpp  $(vpath_src)
vpath %.h    $(vpath_src)

# The X11 base dir on your system
X11BASE=/usr/X11R6
# Add directories with X11 include files here
X11INCS=-I$(X11BASE)/include
# put X11 required libraries and directories here
X11LIBS=-L$(X11BASE)/lib -lX11

SDLDEFS = -D__XWIN__

I_DIRS=-I.. 
#P_DEFS=-DGRAPHICS -DPERIODIC_BOUNDARY

CXXFLAGS = -O3 -Wall -std=gnu++0x $(I_DIRS) $(X11INCS)  $(SDLDEFS) $(P_DEFS)
#CXXFLAGS = -g -Wall -std=gnu++0x $(I_DIRS) $(X11INCS)  $(SDLDEFS) $(P_DEFS)

O = mainMultispeciesMf.o mfSBA.o RWFile.o multiSpecies.o

L = -lm -ltiff

MAIN_TARGET=multiSpeciesSBA

all: $(O)
	g++ -o $(MAIN_TARGET) $(O) $(L)

clean:
	rm $(MAIN_TARGET) $(O) 


# DEPENDENCIES
all:

RWFile.o: RWFile.cpp RWFile.h

mainMultispeciesMf.o: mainMultispeciesMf.cpp 

mfSBA.o: mfSBA.cpp RWFile.cpp RWFile.h

multiSpecies.o: multiSpecies.cpp
