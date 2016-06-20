CXX:=g++
CXXFLAGS:=-O3
EXE:=example.exe
LD_LIB_FLAGS=-L./
LD_INC_FLAGS=-I./


gfa_sort: gfa_sort.o lib
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

gfa_sort.o: gfa_sort.cpp
	$(CXX) $(CXXFLAGS) -c $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

gfa_stats: gfa_stats.cpp lib
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

lib: gfakluge.o
	ar -rs libgfakluge.a $<

test: main.o lib
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge
	./test
	diff q_redundant.gfa q_test.gfa
	cat test.gfa | sort > x.sort; cat test_test.gfa | sort > y.sort; diff x.sort y.sort

main.o: main.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS)

gfakluge.o: gfakluge.cpp gfakluge.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS)

.PHONY: clean

clean:
	$(RM) $(EXE)
	$(RM) *.o
	$(RM) *.a
	$(RM) test
	$(RM) x.sort
	$(RM) y.sort
	$(RM) test_test.gfa
	$(RM) q_test.gfa
	$(RM) gfa_sort
