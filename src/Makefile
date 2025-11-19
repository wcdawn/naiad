# Makefile

CXX = g++
CXXFLAGS = -Wall -Wextra -O3 -std=c++23
#CXXFLAGS = -Wall -Wextra -O0 -g -std=c++23

OBJ = main.o \
	  dense_matrix.o \
	  geometry.o \
	  input.o \
	  output.o \
	  xslibrary.o

EXEC = naiad.x

all : $(EXEC)

$(EXEC) : $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	$(RM) $(OBJ)
	$(RM) $(EXEC)

main.o : input.hpp output.hpp geometry.hpp xslibrary.hpp
input.o : input.hpp geometry.hpp xslibrary.hpp
output.o : output.hpp
geometry.o : geometry.hpp
xslibrary.o : xslibrary.hpp
dense_matrix.o : dense_matrix.hpp
