CXX:=gcc
CXXFLAGS:=-O3 -fopenmp
EXE:=example.exe

$(EXE): main.o gfakluge.o
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS)

gfakluge.o: gfakluge.cpp gfakluge.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS)

.PHONY: clean

clean:
	$(RM) $(EXE)
	$(RM) *.o
