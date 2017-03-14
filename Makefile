CXX?=g++
CXXFLAGS:=-O3 -g -std=c++11
EXE:=example.exe
LD_LIB_FLAGS=-L./src/ -L./
LD_INC_FLAGS=-I./src/ -I./


gfa_sort: src/gfa_sort.cpp lib
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

gfa_stats: src/gfa_stats.cpp lib
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

gfa_verify: src/gfa_verify.cpp lib
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

gfa_diff: src/gfa_diff.cpp lib
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

gfa_merge: src/gfa_merge.cpp lib
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

gfa_ids: src/gfa_ids.cpp lib
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

gfa_spec_convert: src/gfa_spec_convert.o lib
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge
lib: src/gfakluge.o
	ar -rs libgfakluge.a $<

test: src/main.o lib
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge
	./test
	diff q_redundant.gfa q_test.gfa
	cat test.gfa | sort > x.sort; cat test_test.gfa | sort > y.sort; diff x.sort y.sort

src/main.o: src/main.cpp lib
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS)

src/gfakluge.o: src/gfakluge.cpp src/gfakluge.hpp
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
