CXX?=g++
CXXFLAGS:=-O3 -std=c++11 -g
EXE:=example.exe
LD_LIB_FLAGS=-L./src/ -L./
LD_INC_FLAGS=-I./src/ -I./

BIN_DIR:=bin

all: $(BIN_DIR)/gfa_sort $(BIN_DIR)/gfa_subset $(BIN_DIR)/gfa_stats $(BIN_DIR)/gfa_diff $(BIN_DIR)/gfa_merge $(BIN_DIR)/gfa_ids $(BIN_DIR)/gfa_spec_convert $(BIN_DIR)/gfa_extract lib

$(BIN_DIR)/gfa_sort: src/gfa_sort.cpp lib .GFAK_pre-build
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

$(BIN_DIR)/gfa_stats: src/gfa_stats.cpp lib .GFAK_pre-build
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

$(BIN_DIR)/gfa_verify: src/gfa_verify.cpp lib .GFAK_pre-build
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

$(BIN_DIR)/gfa_diff: src/gfa_diff.cpp lib .GFAK_pre-build
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

$(BIN_DIR)/gfa_merge: src/gfa_merge.cpp lib .GFAK_pre-build
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

$(BIN_DIR)/gfa_ids: src/gfa_ids.cpp lib .GFAK_pre-build
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

$(BIN_DIR)/gfa_spec_convert: src/gfa_spec_convert.cpp lib .GFAK_pre-build
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

$(BIN_DIR)/gfa_extract: src/gfa_extract.cpp lib .GFAK_pre-build
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

$(BIN_DIR)/gfa_subset: src/gfa_subset.cpp lib .GFAK_pre-build
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

lib: src/gfakluge.o
	ar -rs libgfakluge.a $<

test: src/main.o lib
	$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge
	./test
	diff data/q_redundant.gfa data/q_test.gfa
	cat data/test.gfa | sort > data/x.sort; cat data/test_test.gfa | sort > data/y.sort; diff data/x.sort data/y.sort

src/main.o: src/main.cpp lib
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS)

src/gfakluge.o: src/gfakluge.cpp src/gfakluge.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS)

.PHONY: clean all


.GFAK_pre-build: 
	mkdir -p bin
	touch .GFAK_pre-build
clean:
	$(RM) $(EXE)
	$(RM) *.o
	$(RM) *.a
	$(RM) test
	$(RM) x.sort
	$(RM) y.sort
	rm -rf $(BIN_DIR)
	$(RM) .GFAK_pre-build
	$(RM) test_test.gfa
	$(RM) q_test.gfa
