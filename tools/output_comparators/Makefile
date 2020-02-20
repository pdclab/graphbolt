#INTT = -DLONG
#INTE = -DEDGELONG

#compilers
$(info ************  Using CILK ************)
PCC = g++-5
PCFLAGS = -std=c++14 -fcilkplus -lcilkrts -g -O3 -DCILK $(INTT) $(INTE)

# dependencies
COMMON = ../../core/common/utils.h ../../core/common/parseCommandLine.h ../../core/common/parallel.h ../../core/common/quickSort.h ../../core/common/blockRadixSort.h ../../core/common/transpose.h ../../core/common/rwlock.h
LOCAL_COMMON = ../common/graphIO.h ../common/fileUtils.h

OUTPUT_COMPARATORS = CFOutputComparator PROutputComparator LPOutputComparator

.PHONY: all clean

all: $(OUTPUT_COMPARATORS)

% : %.C $(COMMON) $(LOCAL_COMMON)
	$(PCC) $(PCFLAGS) -o $@ $<

clean :
	rm -f *.o $(OUTPUT_COMPARATORS)
