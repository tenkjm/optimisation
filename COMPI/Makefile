# To make run
# make dep
# make

# C++ compiler
CPP = c++

# Snowgoose utility folder
SNOWGOOSE = ../snowgoose

# include folders
INCLUDE = $(SNOWGOOSE)

# C++ compiler options
CPPOPTS = --std=c++11 -I$(INCLUDE)

# Libraries to inculde 
LIBS = 

TESTS = testmp.exe testfuncproj.exe

all: testmp.exe testfuncproj.exe

-include deps.inc

clean: 
	rm -f *.exe *.o *.a *~ *.log deps.inc

dep:
	$(CPP) $(CPPOPTS) -MM -c *.cpp >> deps.inc;\
        true
tests:
	@for i in $(TESTS); do if ./$$i > /dev/null; then echo TEST PASSED; continue; else echo TEST FAILED; fi done


.o.exe:
	$(CPP) $(OPTS) -o $@ $< $(LIBS) $(EXTLIBS) $(LDFLAGS)

.cpp.o:
	$(CPP) $(CPPOPTS) -c $<

.c.o:
	$(CC) $(COPTS) -c $<

.SUFFIXES:
.SUFFIXES: .o .a .cpp .c .exe

