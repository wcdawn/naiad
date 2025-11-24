# Makefile

CXX = g++
CXXFLAGS = -Wall -Wextra -O3 -std=c++23
#CXXFLAGS = -Wall -Wextra -O0 -g -std=c++23

OBJ = main.o \
	  dense_matrix.o \
	  exception_handler.o \
	  geometry.o \
	  input.o \
	  output.o \
	  xslibrary.o \
	  diffusion.o \
	  result.o \
	  linalg.o

EXEC = naiad.x

all : $(EXEC)

$(EXEC) : $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	$(RM) $(OBJ)
	$(RM) $(EXEC)

main.o : input.hpp output.hpp geometry.hpp xslibrary.hpp exception_handler.hpp diffusion.hpp result.hpp
input.o : input.hpp geometry.hpp xslibrary.hpp exception_handler.hpp
output.o : output.hpp
geometry.o : geometry.hpp exception_handler.hpp
xslibrary.o : xslibrary.hpp exception_handler.hpp
dense_matrix.o : dense_matrix.hpp
exception_handler.o : exception_handler.hpp
diffusion.o : diffusion.hpp input.hpp geometry.hpp xslibrary.hpp result.hpp exception_handler.hpp linalg.hpp
result.o : result.hpp
linalg.o : linalg.hpp
