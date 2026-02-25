# Makefile

CXX = g++
CFLAGS = -Wall -Wextra -std=c++23 -O3

OBJ = main.o \
	  quadrature1d.o \
	  quadrature_gauss_legendre.o

EXEC = quad.x

.PHONY : default all clean

default : all
all : $(EXEC)

$(EXEC) : $(OBJ)
	$(CXX) $(CFLAGS) -o $@ $^

%.o : %.cpp
	$(CXX) $(CFLAGS) -c $<

main.o : quadrature1d.hpp quadrature_gauss_legendre.hpp
quadrature1d.o : quadrature1d.hpp
quadrature_gauss_legendre.o : quadrature_gauss_legendre.hpp quadrature1d.hpp

clean:
	$(RM) $(OBJ)
	$(RM) $(EXEC)
