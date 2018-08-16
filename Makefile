CPPFLAGS=-Isrc/tinyFA
CXXFLAGS=-std=c++11 -g -O3 -Wall

all: gfak

clean:
	rm -f gfak libgfakluge.a src/*.o

.PHONY: all clean

gfak: src/main.o libgfakluge.a
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

libgfakluge.a: src/gfakluge.o
	$(AR) -rs $@ $^
