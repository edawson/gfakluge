CXX?=g++
CXXFLAGS += -O3 -pipe -fPIC -Wall -std=c++11 -ggdb
PREFIX=/usr/local

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

gfak: $(BUILD_DIR)/main.o src/gfakluge.hpp src/tinyFA/pliib.hpp src/tinyFA/tinyfa.hpp | $(BUILD_DIR) $(BIN_DIR)
	+$(CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS)

$(BUILD_DIR)/main.o: src/main.cpp src/gfakluge.hpp src/tinyFA/pliib.hpp src/tinyFA/tinyfa.hpp | $(BUILD_DIR) $(BIN_DIR)
	+$(CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

install: gfak
	mkdir -p $(DESTDIR)$(PREFIX)/include
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	cp src/gfakluge.hpp $(DESTDIR)$(PREFIX)/include/
	cp src/tinyFA/tinyfa.hpp $(DESTDIR)$(PREFIX)/include/
	cp src/tinyFA/pliib.hpp $(DESTDIR)$(PREFIX)/include/

	cp gfak $(DESTDIR)$(PREFIX)/bin

install-local: gfak
	mkdir -p $(DESTDIR)/include
	mkdir -p $(DESTDIR)/bin
	cp src/gfakluge.hpp $(DESTDIR)/include/
	cp src/tinyFA/tinyfa.hpp $(DESTDIR)/include/
	cp src/tinyFA/pliib.hpp $(DESTDIR)/include/

	cp gfak $(DESTDIR)/bin

check : gfak
	prove test/gfa_test.t


.PHONY: clean all install check


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
