# Makefile for compiling Bowers' Sow & Grow, a C++ project with OpenMP support

# Compiler and compilation flags
CC = g++
CFLAGS = -c -fopenmp -O3
LDFLAGS =
LIBS = -fopenmp -O3

# Project name and executable name
PROJ = bsng
APP = $(PROJ)

# Source files, header files, and object files
SRCS = $(wildcard *.cpp)
HDRS = $(wildcard *.h)
OBJS = $(SRCS:.cpp=.o)

# Targets and their rules
all: $(APP)

$(APP): $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $(APP) $(LIBS)

%.o: %.cpp $(HDRS)
	$(CC) $(CFLAGS) $< -o $@

# Clean up compiled files
clean:
	rm -f *.o $(APP)


#Instructions for Use
#   make        -- compiles your project into program.exe
#   make clean  -- removes compiled item
#   make handin -- creates a project Zip file for hand in
#   All .cpp flles are included.