# Enable this if the input graph has weights. Disable this if you want to use synthetically generated weights.
ifdef WEIGHTED
EDGEDATA = -DEDGEDATA
endif

# Enable this if number of vertices > 2^32 
ifdef LONGVERTEXCOUNT
INTV = -DLONG
endif

INTE = -DEDGELONG

#compilers
$(info ************  Using CILK ************)
PCC = g++
PCFLAGS = -std=c++14 -fcilkplus -lcilkrts -g -O3 -DCILK $(INTV) $(INTE) $(EDGEDATA)
LDFLAGS = -L../lib/mimalloc/out/release -lmimalloc 

# dependencies
COMMON=../core/common/binary_search.h ../core/common/bitsetscheduler.h ../core/common/blockRadixSort.h ../core/common/densebitset.h ../core/common/gettime.h ../core/common/index_map.h ../core/common/matrix.h ../core/common/maybe.h ../core/common/parallel.h ../core/common/parseCommandLine.h ../core/common/transpose.h ../core/common/quickSort.h ../core/common/rwlock.h ../core/common/sequence.h ../core/common/utils.h

GRAPH= ../core/graph/graph.h ../core/graph/graphUtils.h ../core/graph/IO.h ../core/graph/vertex.h ../core/graph/vertexSubset.h

GRAPHBOLT=../core/graphBolt/AdaptiveExecutor.h ../core/graphBolt/GraphBoltEngine_complex.h ../core/graphBolt/GraphBoltEngine.h ../core/graphBolt/GraphBoltEngine_simple.h ../core/graphBolt/ingestor.h ../core/graphBolt/KickStarterEngine.h

OTHERS=../core/main.h

ALL=PageRank LabelPropagation CF COEM SSSP BFS

# make

all: $(ALL)

% : %.C $(COMMON) $(GRAPH) $(GRAPHBOLT) $(OTHERS)
	$(PCC) $(PCFLAGS) $< $(LDFLAGS) -o $@

.PHONY : clean

clean :
	rm -f *.o *.obj $(ALL)

