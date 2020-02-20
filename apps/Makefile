
#INTT = -DLONG
#INTE = -DEDGELONG

#compilers
$(info ************  Using CILK ************)
PCC = g++-5
PCFLAGS = -std=c++14 -fcilkplus -lcilkrts -g -O3 -DCILK $(INTT) $(INTE)
LDFLAGS = -L../lib/mimalloc/out/release -lmimalloc

#

# dependencies
COMMON=../core/common/binary_search.h ../core/common/bitsetscheduler.h ../core/common/blockRadixSort.h ../core/common/densebitset.h ../core/common/gettime.h ../core/common/index_map.h ../core/common/matrix.h ../core/common/maybe.h ../core/common/parallel.h ../core/common/parseCommandLine.h ../core/common/transpose.h ../core/common/quickSort.h ../core/common/rwlock.h ../core/common/sequence.h ../core/common/utils.h

GRAPH= ../core/graph/graph.h ../core/graph/graphUtils.h ../core/graph/IO.h ../core/graph/vertex.h ../core/graph/vertexSubset.h

GRAPHBOLT=../core/graphBolt/AdaptiveExecutor.h ../core/graphBolt/GraphBoltEngine_complex.h ../core/graphBolt/GraphBoltEngine.h ../core/graphBolt/GraphBoltEngine_simple.h ../core/graphBolt/ingestor.h ../core/graphBolt/KickStarterEngine.h

GRAPHBOLT_DEBUG=../core/graphBolt/AdaptiveExecutor.h ../core/graphBolt/debug/GraphBoltEngine_complex_debug.h ../core/graphBolt/debug/GraphBoltEngine_debug.h ../core/graphBolt/debug/GraphBoltEngine_simple_debug.h ../core/graphBolt/ingestor.h ../core/graphBolt/KickStarterEngine.h ../core/graphBolt/ingestor_old.h

OTHERS=../core/main.h

ALL=PageRank LabelPropagation CF SSSP BFS

# make

all: $(ALL)

% : %.C $(COMMON) $(GRAPH) $(GRAPHBOLT) $(OTHERS)
	$(PCC) $(PCFLAGS) $< $(LDFLAGS) -o $@

.PHONY : clean

clean :
	rm -f *.o *.obj $(ALL)

