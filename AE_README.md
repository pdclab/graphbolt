# Artifact Evaluation

Below is a time estimate for different tasks in this readme. The compute time (referred as `compute-seconds`) is based on the provided dataset.

- 1. Overview (5 human-minutes)
- 2. Setting up GraphBolt (10 human-minutes + 10 compute-seconds)
    - 2.1 Requirements (2 human-minutes + 10 compute-seconds if g++/cmake already installed)
    - 2.2 Compiling Applications (2 human-minutes + 10 compute-seconds)
    - 2.3 Dataset Details (5 human-minutes)
        - 2.3.1 Preparing Streaming Datasets for any Graph (Advice to Reviewer: Skip this time-consuming step by directly using provided Wiki-Vote dataset)
- 3. Running Applications (10 human-minutes + 5 compute-seconds)
- 4. Experiments in Paper (20 human-minutes + 15 compute-seconds for an execution from each section)
    - 4.1 Hardware Used for Evaluation (3 human-minutes)
    - 4.2 Evaluating the System (17 human-minutes + 15 compute-seconds for an execution from each section)
        - 4.2.1 Understanding the Output (5 human-minutes)
        - 4.2.2 Varying Batch Sizes (3 human-minutes + 1 compute-second per batch size)
        - 4.2.3 Adaptive Execution (3 human-minutes + 1 compute-second per execution)
        - 4.2.4 Number of Edges Processed (3 human-minutes + 11 compute-seconds per execution)
        - 4.2.5 Sensitivity Experiments (3 human-minutes + 1 compute-second per execution)


## 1. Overview

This readme provides instructions for reproducing the experiments in our paper, `Sparsity-Aware Processing of Streaming Graphs`. 

The sparsity-aware incremental processing technique is implemented in the GraphBolt runtime in order to retain efficiency in presence of sparse computations, thereby pushing the boundary of dependency-driven processing of streaming graphs.

The GraphBolt runtime offers several different processing capabilities (e.g., different modes to ingest streaming graphs), details of which are not relevant for artifact evaluation. Hence, this readme only focuses on the necessary parts to help evaluate the artifact. If you are interested in more details of the GraphBolt and DZiG implementation, refer to the [README.md](README.md) file.

##  2. Setting up GraphBolt

### 2.1 Requirements
- g++ >= 5.3.0 with support for Cilk Plus (Note: gcc-5 and gcc-7 come with cilk support by default.)
- cmake
- [Mimalloc](https://github.com/microsoft/mimalloc) - A fast general purpose memory allocator from Microsoft (version 1.6).
    - Use the helper script `install_mimalloc.sh` to install mimalloc.
      ```bash
      sh install_mimalloc.sh
      ```
    - Update the LD_PRELOAD enviroment variable as specified by install_mimalloc.sh script.

**Important: DZiG requires mimalloc to function correctly and efficiently.**

### 2.2 Compiling Applications

The code for all the applications can be found at `apps/` directory. To compile, run
```bash
cd apps
make -j
```

### 2.3 Dataset Details
  The input graphs used in evaluation are: [UK](http://law.di.unimi.it/webdata/uk-2005/), [TW](http://konect.uni-koblenz.de/networks/twitter), [TT](http://law.di.unimi.it/webdata/twitter-2010/), [FT](https://snap.stanford.edu/data/com-Friendster.html) and [YH](https://webscope.sandbox.yahoo.com/?guccounter=1).
  
[Section 2.3.1](#231-preparing-streaming-datasets-for-any-graph) gives the steps to prepare correct streaming inputs from these datasets. 

  All the above input graphs are very large in size (billions of edges), and converting them to the appropriate format and creating input streams is a time consuming process. To simplify the artifact evaluation process, we have provided the input streams for a smaller graph called [Wiki-Vote](https://snap.stanford.edu/data/wiki-Vote.html) in the `inputs/wiki_vote/`. Hence, you can skip section 2.3.1 and directly start with [section 3](#3-running-applications-1-compute-second-x-5-applications-=-5-compute-seconds) with the Wiki-Vote graph. 

#### 2.3.1 Preparing Streaming Datasets for any Graph

The graph files are obtained in the SNAP format (edge list) format (say, `graph.snap`). In order to create the datasets for our evaluation, do the following:
1. Randomly distribute the lines (excluding any comments) in the graph file to 2 different files as follows:
  ```bash
  # Remove any comments in the graph.snap file
  # Shuffle the lines in the snap file to ensure that the graph is not split based on any order.
  shuf graph.snap -o graph_shuffled.snap 
  # Obtain the line numbers in the snap file and use it to divide it into 2 files
  wc -l graph_shuffled.snap 
  head -n 1000 graph_shuffled.snap > initial_graph.snap
  tail -n 1000+1 graph_shuffled.snap > additions.snap
  ```
2. Convert the snap file to adjacency list format as follows
  ```bash
  cd tools/converters
  make SNAPtoAdjConverter
  ./SNAPtoAdjConverter initial_graph.snap initial_graph.adj
  # for undirected (symmetric) graphs, use the -s flag
  ./SNAPtoAdjConverter -s initial_graph.snap initial_graph.adj.un
  ```
  More details in section 2.4 of [README.md](README.md#24-graph-input-and-stream-input-format).

3. Creating the streaming input
  * The `additions.snap` file is used to create edge additions stream file, `additions.stream` as follows:
    ```bash
    sed -e 's/^/a\t/' additions.snap > additions.stream
    ```
  * The `initial_graph.snap` file is used to create edge deletions stream file, `deletions.stream`. This is one of the steps to ensure that all the deletion operations actually result in edge deletions. In addition to this, you also need to use the command-line parameters `-fixedBatchSize -enforceEdgeValidity` while running the application as described in the next [section](#3-running-applications-1-compute-seconds-x-5-applications-=-5-compute-seconds):
    ```bash
    sed -e 's/^/d\t/' initial_graph.snap > deletions.stream
    ```
  * Then to combine the additions and deletions file into the same file, the following command can be used:
    ```bash
    paste -d "\n" additions.stream deletions.stream > update.stream
    ```
4. CF and COEM applications work on bipartite input graphs. In order to mock bipartite graphs using non-bipartite graphs, we create a partitions file which contains all the vertices belonging to a single partition. We do this by randomly assigning each vertex in the graph to one of the partitions.
    ```bash
    cd tools/generators
    # For an input graph with 8298 vertices, do the following.
    python CreatePartitions.py --n=8298 --outputFile="Partition1_wiki"
    ```
5. Label Propagation and COEM use an initial set of seed vertices as the frontier. We randomly select 10% of the vertices as the seeds for each input graph. The seeds file can be generated as follows:
    ```bash
    cd tools/generators
    # For an input graph with 8298 vertices, do the following.
    python CreateSeeds.py --n=8298 --outputFile="LabelSeeds_wiki" --seedsPercent=10
    ```


### 3. Running Applications (1 compute-second x 5 applications = 5 compute-seconds)

The command for running each application with the `wiki_vote` graph is provided below:
```bash
# Running PageRank
./PageRank -maxIters 10 -fixedBatchSize -nEdges 1000 -streamPath ../inputs/wiki_vote/wiki_valid.stream -numberOfUpdateBatches 2 ../inputs/wiki_vote/wiki_vote_initial.adj

# Running Label Propagation
./LabelPropagation -maxIters 10 -fixedBatchSize -seedsFile ../inputs/wiki_vote/LabelSeeds_wiki -features 2 -nEdges 1000 -streamPath ../inputs/wiki_vote/wiki_valid.stream -numberOfUpdateBatches 2 ../inputs/wiki_vote/wiki_vote_initial.adj

# Running COEM
./COEM -s -maxIters 10 -fixedBatchSize -seedsFile ../inputs/wiki_vote/LabelSeeds_wiki -partitionsFile ../inputs/wiki_vote/Partition1_wiki -nEdges 1000 -streamPath ../inputs/wiki_vote/wiki_valid_sym_partitioned.stream -numberOfUpdateBatches 2 ../inputs/wiki_vote/wiki_vote_initial.adj.un

# Running CF
./CF -s -maxIters 10 -fixedBatchSize -partitionsFile ../inputs/wiki_vote/Partition1_wiki -modVal 0.022000001 -numberOfFactors 2 -nEdges 1000 -streamPath ../inputs/wiki_vote/wiki_valid_sym_partitioned.stream -numberOfUpdateBatches 2 ../inputs/wiki_vote/wiki_vote_initial.adj.un

# Running SSSP
./SSSP -source 2565 -weightCap 5 -fixedBatchSize -nEdges 1000 -streamPath ../inputs/wiki_vote/wiki_valid.stream -numberOfUpdateBatches 2 ../inputs/wiki_vote/wiki_vote_initial.adj
```

**IMPORTANT: Ensure that `-fixedBatchSize` flag is present.** This ensures that the specified number of edge operations `-nEdges` are correctly performed. In addition to this, `-enforceEdgeValidity` flag validates all edge operations before counting them towards the batch size requirement (for example, deleting edges not present in the graph is considered as an invalid edge operation). For further details regarding ingestor flags, refer [Stream Ingestors](README.md#5-Stream-Ingestor) section in the [README.md](README.md) file.

For our experiments, we first create the validated additions and deletions files for each input graph. And then, combine these files to get the validated update stream. This ensures that for any given batch size, there are equal number of additions and deletions. The datasets for the provided `Wiki-Vote` graph are already validated, and hence `-enforceEdgeValidity` is not needed when running the applications. 

##  4. Experiments in Paper

### 4.1 Hardware Used for Evaluation
All the experiments in the paper are carried out in the following 2 machines:

`System 1` - All the experiments (except those on Yahoo dataset) are carried out in the 32 core machine running a 2GHz processor with 231 GB RAM. Details of this machine can be found [here](system1_environment.log). This information is obtained by running the `collect_environment.sh` script from [here](https://github.com/SC-Tech-Program/Author-Kit).

`System 2` - Experiments on the Yahoo dataset are carried out in the r5.24xlarge machine on Amazon EC2. This is a 96 core machine with 748GB RAM running at 2.5GHz. Details of this machine can be found [here](system2_environment.log). 


### 4.2 Evaluating the System

#### 4.2.1 Understanding the Output
Below is a sample output of running PageRank with `wiki_vote` graph: 
```bash
Graph created
Initializing engine ....
Number of batches: 1
Creating dependency structure ....
Initializing dependency structure ....
Finished initializing engine

   Iteration,        Time
           1,    0.000475
           2,    0.000417
           3,    0.000285
           4,    0.000284
           5,    0.000234
           6,    0.000238
           7,    0.000223
           8,    0.000218
           9,    0.000218
          10,    0.000189
Initial graph processing : 0.003682
Number of iterations : 10

Opening Stream: Waiting for writer to open...
Stream opened
Current_batch: 1
Batch Size: 10000
Reading Time : 0.007327
Edge deletion time : 0.000898
Edge addition time : 0.000205
Edge Additions in batch: 5000
Edge Deletions in batch: 5000

   Iteration,        Time
           1,    0.000503
           2,    0.000372
           3,    0.000278
           4,    0.000312
           5,    0.000301
           6,    0.000281
           7,    0.000263
           8,    0.000273
           9,    0.000253
          10,    0.000225
Finished batch : 0.003281
Number of iterations : 10

Hit Max Batch Size
```
The `Finished batch` line gives the execution time of DZiG in seconds to process each batch of edge updates. In the above output sample, DZiG took 0.003281 seconds. The execution times of `DZiG` throughout the paper (including Figure 10, 11, 13, 15, 17 and 18) are obtained from the `Finished batch` line.

The time taken in seconds to apply the edge additions and edge deletions is printed as `Edge addition time` and `Edge deletion time`. In the above output sample, it took 0.000898 seconds to apply edge deletions and 0.000205 seconds to apply edge additions. Graph mutation times reported in Figure 18 are obtained by summing the addition and deletion times.

#### 4.2.2 Varying Batch Sizes 

Throughout our evaluation, we presented results for different applications running with different input batch sizes. This is achieved by varying the `-numberOfUpdateBatches` command-line parameter.

#### 4.2.3 Adaptive Execution

In order to enable adaptive execution, use the command-line parameter `-ae`. For example, 
```bash
./PageRank -ae -maxIters 10 -fixedBatchSize -nEdges 10000 -streamPath ../inputs/wiki_vote/wiki_valid.stream -numberOfUpdateBatches 2 ../inputs/wiki_vote/wiki_vote_initial.adj
```
The execution times for `DZiG-AE` is used in Figure 15.

#### 4.2.4 Number of Edges Processed

Figures 12, 14 and 16 show the amount of work done in terms of the number of edges processed by DZiG. This is obtained by performing the following steps.

Compile the application with `WORKFLAG=1` as follows:
```bash
make clean 
make WORKFLAG=1 PageRank
```

And then, run the application:
```bash
./PageRank -maxIters 10 -fixedBatchSize -nEdges 10000 -streamPath ../inputs/wiki_vote/wiki_valid.stream -numberOfUpdateBatches 2 ../inputs/wiki_vote/wiki_vote_initial.adj
```

The output will contain a line with `Edges Processed`, as shown below. 

```bash
Edges Processed : 357552
```

#### 4.2.5 Sensitivity Experiments
In order to vary the epsilon value as done in Figure 18, use the command-line parameter `-epsilon` as follows:
```bash
./PageRank -epsilon 0.001 -maxIters 10 -fixedBatchSize -nEdges 10000 -streamPath ../inputs/wiki_vote/wiki_valid.stream -numberOfUpdateBatches 2 ../inputs/wiki_vote/wiki_vote_initial.adj
```

