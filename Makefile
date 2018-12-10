CXX?=g++
CXXFLAGS:=-O3 -pipe -fPIC -march=native -mtune=native -Wall -std=c++11 -g -ggdb
# We want to pass -Wa,-q to GCC use the Clang assembler, but Apple Clang can't take that
# So we do an environment variable instead
export AS_INTEGRATED_ASSEMBLER=1

ifeq ($(shell if [ -d /opt/local/include/libomp ];then echo 1;else echo 0;fi), 1)
    # On OS X with Apple Clang, <omp.h>, which our tinyFA dependency needs, isn't always on the include path
    # Pick it up from Macports if it is there.
    # Homebrew ought to put it where the compiler can find it.
    CXXFLAGS += -I/opt/local/include/libomp
endif


BIN_DIR:=bin
BUILD_DIR:=build

LD_LIB_FLAGS=-L./src/ -L./
LD_INC_FLAGS=-I./src/ -I./ -I./src/tinyFA -I./$(BUILD_DIR)

gfak: $(BUILD_DIR)/main.o libgfakluge.a src/tinyFA/pliib.hpp src/tinyFA/tinyfa.hpp | $(BUILD_DIR) $(BIN_DIR)
	+$(CXX) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

$(BUILD_DIR)/main.o: src/main.cpp src/gfakluge.hpp src/tinyFA/pliib.hpp src/tinyFA/tinyfa.hpp | $(BUILD_DIR) $(BIN_DIR)
	+$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS) -lgfakluge

libgfakluge.a: $(BUILD_DIR)/gfakluge.o src/tinyFA/pliib.hpp src/tinyFA/tinyfa.hpp | $(BUILD_DIR) $(BIN_DIR)
	ar -rs $@ $<

$(BUILD_DIR)/gfakluge.o: src/gfakluge.cpp src/gfakluge.hpp src/tinyFA/pliib.hpp src/tinyFA/tinyfa.hpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)


.PHONY: clean all


clean:
	$(RM) gfak
	$(RM) src/*.o
	$(RM) *.a
	$(RM) x.sort
	$(RM) y.sort
	rm -rf $(BIN_DIR)
	rm -rf $(BUILD_DIR)
	$(RM) test_test.gfa
	$(RM) q_test.gfa
