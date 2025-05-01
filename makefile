# Makefile for the sequential IST generator (sequential.cpp)
# - builds 'sequential' with -O2 and gprof support (-pg)
# - 'make clean' removes everything except parents.csv
# - 'make run N=8' 
# - 'make gprof N=8' will run the binary and generate a gprof report

CXX       := g++
CXXFLAGS  := -O2 -g 
LDFLAGS   := -pg

SRC       := sequential.cpp
TARGET    := sequential
CSV       := Sequential_parents.csv


.PHONY: all clean run gprof

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

run: all
	@echo "Running $(TARGET) with N=$(N)…"
	@./$(TARGET) $(N)

gprof: all
	@echo "Profiling run with N=$(N)…"
	@./$(TARGET) $(N)
	@echo "Generating gprof report to profile.txt"
	@gprof $(TARGET) gmon.out > profile.txt

clean:
	@echo "Cleaning binaries and profiles (keeping $(CSV))"
	-@rm -f $(TARGET) *.o gmon.out profile.txt
