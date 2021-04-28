# Detailed Instructions for Converting Datasets Used in Artifact Evaluation

## UK Dataset from WebData

1. Download the dataset (we require both the .graph and .properties files):
	```bash
	wget http://data.law.di.unimi.it/webdata/uk-2005/uk-2005.graph
	wget http://data.law.di.unimi.it/webdata/uk-2005/uk-2005.properties
	```
2. Clone this repository to get the `.graph` converter
	```bash
	git clone https://github.com/Mogami95/graph-xll.git
	```
3. Enter the directory of `graph-xll`
	```bash
	cd graph-xll
	```
4. Run converter command. The java command takes as input the base dataset (<DATASET_FOLDER_PATH> is the folder path where the dataset is placed). 
	```bash
	java -cp "lib/*:bin" BV2Ascii <DATASET_FOLDER_PATH>/uk-2005 > <DATASET_FOLDER_PATH>/uk-2005.snap
	```
5. Follow the instructions from [Section 2.3.1 of the AE_README](AE_README.md#231-preparing-streaming-datasets-for-any-graph)



## TW Dataset from Konect

1. Download the dataset:
	```bash
	wget http://konect.cc/files/download.tsv.twitter.tar.bz2
	```
2. Unzip the dataset
	```bash
	tar -xvjf download.tsv.twitter.tar.bz2
	```
3. Remove any header/comments from the file (comments can be either % or #)
	```bash
	sed '/^[#%]/d' out.twitter > out.twitter.snap
	```
4. Follow the instructions from [Section 2.3.1 of the AE_README](AE_README.md#231-preparing-streaming-datasets-for-any-graph)


## TT Dataset from Konect

1. Download the dataset:
	```bash
	wget http://konect.cc/files/download.tsv.twitter_mpi.tar.bz2
	```
2. Unzip the dataset
	```bash
	tar -xvjf download.tsv.twitter_mpi.tar.bz2
	```
3. Remove any header/comments from the file (comments can be either % or #)
	```bash
	sed '/^[#%]/d' out.twitter_mpi > out.twitter_mpi.snap
	```
4. Follow the instructions from [Section 2.3.1 of the AE_README](AE_README.md#231-preparing-streaming-datasets-for-any-graph)



## FT Dataset from Konect

1. Download the dataset:
	```bash
	wget http://konect.cc/files/download.tsv.friendster.tar.bz2
	```
2. Unzip the dataset
	```bash
	tar -xvjf download.tsv.friendster.tar.bz2
	```
3. Remove any header/comments from the file (comments can be either % or #)
	```bash
	sed '/^[#%]/d' out.friendster > out.friendster.snap
	```
4. Follow the instructions from [Section 2.3.1 of the AE_README](AE_README.md#231-preparing-streaming-datasets-for-any-graph)


## YH Dataset

1. Request the `G2 - Yahoo! AltaVista Web Page Hyperlink Connectivity Graph, circa 2002 (multi part) (Hosted on AWS)` from yahoo webscope. It may take a couple of days for the request to be approved.
2. After the request is approved follow the instructions in the approval e-mail to download the dataset.
3. In the artifacts branch of our repo, update the repo to get a new converter
    ```bash
	git pull
	```
4. From the main directory of the repo, enter the converter directory
    ```bash
	cd tools/converters/
	```
5. Create the binary for the converter
    ```bash
	make adjacencyToSNAP
	```
6. The dataset is split into multiple parts. Convert the graph data files (ydata-yaltavista-webmap-v1_0_links-*.txt) into the snap format. In below commands, <DATASET_FOLDER_PATH> is the folder path where the dataset is placed.
    ```bash
	./adjacencyToSNAP --adjacencylistfile=<DATASET_FOLDER_PATH>/ydata-yaltavista-webmap-v1_0_links-1.txt --undirected=0 --header=0
	./adjacencyToSNAP --adjacencylistfile=<DATASET_FOLDER_PATH>/ydata-yaltavista-webmap-v1_0_links-2.txt --undirected=0 --header=0
	```
7. Concat the two files to form a single snap dataset.
    ```bash
	cat <DATASET_FOLDER_PATH>/ydata-yaltavista-webmap-v1_0_links-1.txt.snap <DATASET_FOLDER_PATH>/ydata-yaltavista-webmap-v1_0_links-2.txt.snap > <DATASET_FOLDER_PATH>/ydata-yaltavista-webmap-v1_0_links-fulldataset.snap
	```
8. Follow the instructions from [Section 2.3.1 of the AE_README](AE_README.md#231-preparing-streaming-datasets-for-any-graph)
