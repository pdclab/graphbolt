<img src="https://user-images.githubusercontent.com/8582843/60870031-e6e69f80-a1e4-11e9-8d44-e8472355230a.png" alt="GraphBolt" width="200">

**GraphBolt: Dependency-Driven Synchronous Processing of Streaming Graphs**

## 1. What is it?

GraphBolt is an efficient streaming graph processing system that provides Bulk Synchronous Parallel (BSP) guarantees. GraphBolt performs dependency-driven incremental processing which quickly reacts to graph changes, and provides low latency & high throughput processing. [[Read more]](https://www.cs.sfu.ca/~keval/contents/papers/graphbolt-eurosys19.pdf)

For asynchronous algorithms, GraphBolt incorporates KickStarter's light-weight dependency tracking and trimming strategy. [[Read more]](https://www.cs.sfu.ca/~keval/contents/papers/kickstarter-asplos17.pdf)

##  2. Getting Started

### 2.1 Core Organization

The `core/graphBolt/` folder contains the [GraphBolt Engine](#3-graphbolt-engine), the [KickStarter Engine](#4-kickstarter-engine), and our [Stream Ingestor](#5-stream-ingestor) module. The application/benchmark codes (e.g., PageRank, SSSP, etc.) can be found in the `apps/` directory. Useful helper files for generating the stream of changes (`tools/generators/streamGenerator.C`), creating the graph inputs in the correct format (`tools/converters/SNAPtoAdjConverter.C` - from ligra's codebase), comparing the output of the algorithms (`tools/output_comparators/`) are also provided.

### 2.2 Requirements
- g++ >= 5.3.0 with support for Cilk Plus.
- [Mimalloc](https://github.com/microsoft/mimalloc) - A fast general purpose memory allocator from Microsoft (version >= 1.6).
    - Use the helper script `install_mimalloc.sh` to install mimalloc.
    - Update the LD_PRELOAD enviroment variable as specified by install_mimalloc.sh script.

**Important: GraphBolt requires mimalloc to function correctly and efficiently.**

### 2.3 Compiling and Running the Application

Compilation is done from within `apps` directory. To compile, run
```bash
$   cd apps
$   make -j
```
 The executable takes the following command-line parameters:
 - `-s` : Optional parameter to indicate a symmetric (undirected) graph is used. 
 - `-streamPath` : Path to the input stream file or pipe (More information on the input format can be found in [Section 2.4](#24-graph-input-and-stream-input-format)).
 - `-numberOfUpdateBatches` : Optional parameter to specify the number of edge updates to be made. Default is 1.
 - `-nEdges` : Number of edge operations to be processed in a given update batch.
 - `-outputFile` : Optional parameter to print the output of a given algorithms.
 - Input graph file path (More information on the input format can be found in [Section 2.4](#24-graph-input-and-stream-input-format)).

For example, 
```bash
$   # Ensure that LD_PRELOAD is set as specified by the install_mimalloc.sh
$   ./PageRank -numberOfUpdateBatches 2 -nEdges 1000 -streamPath ../inputs/testInputFile -outputFile /tmp/output/pr_output ../inputs/testGraph.adj
$   ./SSSP -source 0 -numberOfUpdateBatches 1 -nEdges 1000 -streamPath GRAPHBOLT_HOME/myFifo -outputFile /tmp/output/sssp_output ../inputs/testGraph.adj
```
Other additional parameters may be required depending on the algorithm. Additional configurations for the graph ingestor and the graph can be found in [Section 5](#5-stream-ingestor).

### 2.4 Graph Input and Stream Input Format

The initial input graph should be in the [adjacency graph format](http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html). 
For example, the SNAP format (edgelist) and the adjacency graph format for a sample graph are shown below.
SNAP format:
```txt
0 1
0 2
2 0
2 1
```
 Adjacency Graph format:
```bash
AdjacencyGraph
3
4
0
2
2
1
2
0
1
```
You can use `tools/converters/SNAPtoAdjConverter` to convert an input graph in Edgelist format (SNAP format) to the adjacency graph format, as follows:
```bash
$ ./SNAPtoAdjConverter inputGraph.snap inputGraph.adj
# for undirected (symmetric) graphs, use the -s flag
$ ./SNAPtoAdjConverter -s inputGraph.snap inputGraphUndirected.adj 
```
The streaming input file should have the edge operation (addition/deletion) on a separate line. The edge operation should be of the format, `[d/a] source destination` where `d` indicates edge deletion and `a` indicates edge addition. Example streaming input file:
```bash
a 1 2
d 2 3
a 4 5
...
```

Edge operations can be streamed through a pipe using `tools/generators/streamGenerator.C`. It takes in the following command-line parameters:
- `-edgeOperationsFile` : Input file containing the edge operations in the format mentioned above.
- `-outputPipe` : Path of the output pipe where the edges are streamed to.

```bash
$   cd tools/generators
$   make streamGenerator
$   ./streamGenerator -edgeOperationsFile edgeOperationsFilepath -outputPipe outputPipePath
```
More details regarding the ingestor can be found in [Section 5](#5-stream-ingestor).

## 3. GraphBolt Engine

The GraphBolt engine provides Bulk Synchronous Parallel (BSP) guarantees while incrementally processing streaming graphs.

### 3.1 Creating Applications using the GraphBolt Engine

A key design decision of the GraphBolt framework is to ensure that the application code remains oblivious to GraphBolt's internal subtleties while still providing fast performance.

So, the application code only needs to express its computation using the following functions. More details regarding these functions can be found in the inline comments of `GraphBoltEngine.h`.

#### AggregateValue and VertexValue initialization:
- initializeAggregationValue()
- initializeVertexValue()
- aggregationValueIdentity()
- vertexValueIdentity()

GraphBolt stores information for each vertex in the form of aggregation values. So, first, the user should identify the aggregation value and the vertex value for the algorithm. For example in PageRank, the vertex value is its pagerank (`PR`) and the aggregation value is the sum of `(PR[u]/out_Degree[u])` values from all its inNeighbors. 

#### Activate vertex / Compute vertex for a given iteration:
- forceActivateVertexForIteration()
- forceComputeVertexForIteration()
- shouldUseDelta()

In iterative graph algorithms, at a given iteration `i`, a set of vertices will push some value to their outNeighbors. These are the active vertices for that iteration. The outNeighbors which receive these values will then compute their updated values. The following functions are provided to force a vertex to be either active/compute at a given iteration. For example, in Label Propagation, all the vertices should compute their values at each iteration irrespective of whether they receive any new changes from their inNeighbors at that iteration (refer `apps/LabelPropagation.C`).

#### Add to or remove from aggregation:
- addToAggregation()
- addToAggregationAtomic()
- removeFromAggregation()
- removeFromAggregationAtomic()

These are the functions used to add a value to or remove some value from the aggregation value. For sum, it is simply adding and subtracting the values from the aggregation value passed. Note that `addToAggregationAtomic()` and `removeFromAggregationAtomic()` will be called by multiple threads on the same aggregation value. So, the update should be performed atomically using CAS.

#### Edge functions:
- sourceChangeInContribution()
- edgeFunction()
- edgeFunctionDelta()

The edge operation is split into 3 phases:
1. Determine the source contribution - The computations for a given vertex which are dependent only on the source values are performed here. For example, in PageRank, a vertex `u` adds the value `PR[u]/out_degree[u]` to the aggregation value of all its outNeighbors. Since this computation of `PR[u]/out_degree[u]` is common for processing all the outEdges of `u`, we can compute this value (contribution of the source vertex) only once and perform the addition for all outEdges.
2. Transform the contribution depending on the edge data - In this step, the source vertex contribution is transformed by the edge property. For example in weighted page rank, the contribution will be multiplied by the edge weight.
3. Aggregating the contribution to the aggregation value using `addToAggregationAtomic()`.

Note that these functions do not require CAS or locks. In the case of complex aggregations, an additional `edgeFunctionDelta()` has to be defined. Refer the `apps/GraphBoltEngine.h`, `apps/GraphBoltEngine_complex.h` for further details of these functions.

#### Vertex compute function and determine end of computation:
- computeFunction()
- isChanged()

Given an aggregation value, `computeFunction()` computes the vertex value corresponding to this aggregation value.
In order to detemine the convergence condition, the `isChanged()` is used to determine whether the value of vertex has significantly changed compared to its previous value. 
Both these functions do not require CAS or locks as they will be invoked in a vertex parallel manner.

#### Determine how an edge update affects the source / destination:
- hasSourceChangedByUpdate()
- hasDestinationChangedByUpdate()

Should return true if the source or destination vertex of an edge operation becomes active in the first iteration. For example, in PageRank, if the out_degree of a vertex changes, then it will be active in the first iteration.

#### Compute function
- compute()

This is the starting point of the application. The GraphBolt engine is initialized here with the required configurations and started.

In addition to these functions, the algorithm also needs to define an `Info` class which contains all the global variable/constants required for that application. It should implement the following functions:
- copy()
- processUpdates()
- cleanup()


## 4. KickStarter Engine

The KickStarter engine is used for streaming path-based/monotonic graph algorithms like SSSP, BFS etc.

### 4.1 Creating Applications using KickStarter Engine

Similar to the GraphBolt engine, the KickStarter engine also provides functions to express the algorithm.
#### VertexValue initialization:
- initializeVertexValue()

#### Activate vertex / Compute vertex for first iteration:
- frontierVertex()

Unlike the GraphBolt engine, we only need to know which of the vertices are active in the first iteration.

#### Edge functions:
- edgeFunction()

For an edge (u, v), compute v's value based on u's value. Return false if the value from u should not be use to update the value of v. Return true otherwise.

### ShouldPropagate:
- shouldPropagate()

The shouldPropagate condition to determine whether the monotonicity of the vertex holds given 2 values depending on the algorithm.

#### Compute function
- compute()

The starting point of the application. The KickStarter engine is initialized here with the required configurations and started.


## 5. Stream Ingestor

The stream ingestor FIFO is specified by `-streamPath`. Edge operations can be written to this FIFO. `-nEdges` specifies the maximum number of edge operations that can be passed to the GraphBolt engine in a single batch. The GraphBolt engine will continue to receive batches of edges from the stream ingestor until either the stream is closed (when there are no more writers to the FIFO) or when `-numberOfBatches` has been exceeded. If the writing end of the FIFO is not opened, the GraphBolt engine (which is the reading end) will block and wait until it is opened. 

There are a few optional flags that can affect the behaviour and determine the validity of the edge operations  passed to the command line parameter `-streamPath`:

- `-fixedBatchSize`: Optional flag to ensure that the batch size is strictly adhered to. If the FIFO does not contain enough edges, the ingestor will block until it has received enough edges specified by `-nEdges` or until the stream is closed. 
- `-enforceEdgeValidity`: Optional flag to ensure that all edge operations in the batch are valid. For example, an edge deletion operation is valid only if the edge to be deleted is present in the graph. In the case of a `simple graph` (explained below), an edge addition operation is valid only if that edge does not currently exist in the graph. Invalid edges are discarded and are not included while counting the number of edges in a batch.
- `-simple`: Optional flag used to ensure that the input graph remains a simple graph (ie. no duplicate edges). The input graph is checked to remove all duplicate edges. Duplicate edges are not allowed within a batch and edge additions are checked to ensure that the edge to be added does not yet exist within the graph.
- `-debug`: Optional flag to print the edges that were determined to be invalid.

## 6. Acknowledgements
Some utility functions from [Ligra](https://github.com/jshun/ligra) and [Problem Based Benchmark Suite](http://www.cs.cmu.edu/~pbbs/index.html) are used as part of this project. We are thankful to them for releasing their source code.

## 7. Resources
Mugilan Mariappan and Keval Vora. [GraphBolt: Dependency-Driven Synchronous Processing of Streaming Graphs](https://dl.acm.org/citation.cfm?id=3303974). European Conference on Computer Systems (**EuroSys'19**). Dresden, Germany, March 2019.

Keval Vora, Rajiv Gupta and Guoqing Xu  [KickStarter: Fast and Accurate Computations on Streaming Graphs via Trimmed Approximations](https://dl.acm.org/citation.cfm?id=3093336.3037748). Architectural Support for Programming Languages and Operating Systems (**ASPLOS'17**). Xi'an, China, April 2017.


To cite, you can use the following BibTeX entries:

```
@inproceedings{Mariappan:2019:GDS:3302424.3303974,
 author = {Mariappan, Mugilan and Vora, Keval},
 title = {GraphBolt: Dependency-Driven Synchronous Processing of Streaming Graphs},
 booktitle = {Proceedings of the Fourteenth EuroSys Conference 2019},
 series = {EuroSys '19},
 year = {2019},
 isbn = {978-1-4503-6281-8},
 location = {Dresden, Germany},
 pages = {25:1--25:16},
 articleno = {25},
 numpages = {16},
 url = {http://doi.acm.org/10.1145/3302424.3303974},
 doi = {10.1145/3302424.3303974},
 acmid = {3303974},
 publisher = {ACM},
 address = {New York, NY, USA},
 keywords = {Incremental Processing, Streaming Graphs},
} 

@inproceedings{Vora:2017:KFA:3037697.3037748,
 author = {Vora, Keval and Gupta, Rajiv and Xu, Guoqing},
 title = {KickStarter: Fast and Accurate Computations on Streaming Graphs via Trimmed Approximations},
 booktitle = {Proceedings of the Twenty-Second International Conference on Architectural Support for Programming Languages and Operating Systems},
 series = {ASPLOS '17},
 year = {2017},
 isbn = {978-1-4503-4465-4},
 location = {Xi'an, China},
 pages = {237--251},
 numpages = {15},
 url = {http://doi.acm.org/10.1145/3037697.3037748},
 doi = {10.1145/3037697.3037748},
 acmid = {3037748},
 publisher = {ACM},
 address = {New York, NY, USA},
 keywords = {graph processing, streaming graphs, value dependence},
} 
```