CXX						= g++
CXXFLAGS			= -O3 -g -fomit-frame-pointer  -std=c++14 -Wwrite-strings -pthread
LDFLAGS				= 
LDLIBS				= 
ARG	= log/compres_seq_memory_withoutrawFBWT.dat

# targets
.PHONY: all
all:fbwtmem_seq # fbwtmem_bi 

fbwtmem_seq: sais/sais.o makeFBWT_seq.o fbwtmem_seq.o #fbwtmem_seq.o
	$(CXX) sais/sais.o fbwtmem_seq.o  makeFBWT_seq.o -o fbwtmem_seq  -pthread
fbwtmem_bi: sais/sais.o fbwtmem_bi.o makeFBWT_bi.o
	$(CXX) sais/sais.o fbwtmem_bi.o  makeFBWT_bi.o -o fbwtmem_bi  -pthread

.cpp.o:
	g++ $(CXXFLAGS) -Wall -c $<

distclean: clean
clean:
	$(RM) *.o fbwtmem_seq fbwtmem_bi *~

# dependencies
fwtmem_seq.o: Makefile

# cat $(ARG) | grep memLen | awk '{if(NR % 2 == 0){print $$0}}' | awk '{print $$10,$$12,$$4,$$2,$$18}'
seqLog:
	cat $(ARG) | grep memLen | awk '{if(NR % 2 == 1){print $$0}}' | awk '{print $$10,$$12,$$4,$$18,$$2}'
biLog:
	cat $(ARG) | grep memLen | awk '{if(NR % 2 == 0){print $$0}}' | awk '{print $$10,$$12,$$4,$$2,$$18}'
seqLogOnly:
	bash -c "paste <(cat $(ARG) | grep memLen)  <(cat $(ARG) | grep \"#memory\") | awk '{print \$$10,\$$12,\$$4,\$$18,\$$2,\$$6,\$$14}'"
#	bash -c "paste <(cat $(ARG) | grep memLen)  <(cat $(ARG) | grep \"#memory\") | awk '{print \$$10,\$$12,\$$4,\$$22,\$$2,\$$6,\$$14}'"

# for f in `find ../data/ | grep -e fa -e fna | grep -v Dro`; do a=`basename $f` && b=${a%\.fa} && c=${b%\.fna} && mkdir -p ../index/$c; done
