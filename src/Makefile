# Makefile

CXX = g++
CXXFLAGS = -Wall -Wextra -O3 -std=c++23

OBJ = main.o \
	  input.o \
	  output.o

EXEC = naiad.x

all : $(EXEC)

$(EXEC) : $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

clean:
	$(RM) $(OBJ)
	$(RM) $(EXEC)

main.o : input.hpp output.hpp
input.o : input.hpp
output.o : output.hpp
