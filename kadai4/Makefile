#------------------------------------------------------
# Makefile for C/C++ Program
#------------------------------------------------------
# Target: a.out
# Author: H. Watanabe
#------------------------------------------------------

TARGET=a.out

#------------------------------------------------------
# Default Parameters
#------------------------------------------------------

CC=g++
OPT=
INC=
LIB=

#------------------------------------------------------
# Compile Option
#------------------------------------------------------

-include makefile.opt

#------------------------------------------------------
# Definition
#------------------------------------------------------

.SUFFIXES:.cc .c .o .h

#---
# Source Files
#---

SRC=$(shell ls *.cpp)
HED=$(shell ls *.h)
OBJ=$(SRC:.cpp=.o)

#------------------------------------------------------
# rules
#------------------------------------------------------

all: $(TARGET)
$(TARGET): $(OBJ)
	$(CC) $(OPT) -o $(TARGET) $(OBJ) $(LIB)
	rm *.o

.c.o:
	$(CC) $(OPT) -c $< $(INC)
.cc.o:
	$(CC) $(OPT) -c $< $(INC)

dep:
	g++ -MM -MG $(SRC) >makefile.depend

clean:
	rm -f $(TARGET) $(TARGET).exe
	rm -f *.o *.obj
	rm -f *~ *.~*

tar:
	tar cvzf $(TARGET).tar.gz $(SRC) $(HED) makefile

#--------------------------------------------------
-include makefile.depend