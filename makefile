# Makefile

IDIR = include
CC = g++
CFLAGS = -I$(IDIR) -std=c++14 -Wall  -g -O3 -fopenmp -lboost_program_options

SDIR = src
ODIR = obj

_SOURCES = driver.cpp graph.cpp node.cpp edge.cpp summ.cpp utils.cpp
SOURCES = $(patsubst %,$(SDIR)/%,$(_SOURCES))

_OBJECTS = $(_SOURCES:.cpp=.o)
OBJECTS = $(patsubst %,$(ODIR)/%,$(_OBJECTS))

_DEPS = graph.h node.h edge.h summ_opts.h summ.h utils.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

EXE = summ

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXE): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS)
	rm -f $(ODIR)/*.o
