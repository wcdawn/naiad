# Makefile

CXX = g++
CFLAGS = -Wall -Wextra -std=c++23 -O3

OBJ = quadrature1d.o main.o

EXEC = quad.x

.PHONY : default all clean

default : all
all : $(EXEC)

$(EXEC) : $(OBJ)
	$(CXX) $(CFLAGS) -o $@ $^

%.o : %.cpp
	$(CXX) $(CFLAGS) -c $<

main.o : quadrature1d.hpp
quadrature1d.o : quadrature1d.hpp ../gauss_legendre.txt

clean:
	$(RM) $(OBJ)
	$(RM) $(EXEC)
