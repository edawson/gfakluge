CXX:=g++
CXXFLAGS:= -O3 -fPIC -std=c++11 -pipe -D_FILE_OFFSET_BITS=64
LD_INC_FLAGS:= -I./

getseq: examples/getseq.cpp pliib.hpp tinyfa.hpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_INC_FLAGS)

indexseq: examples/index.cpp pliib.hpp tinyfa.hpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LD_INC_FLAGS)
