__INPUT FILES FORMAT__

__EMB_min__ loads the graph in __edge list__ format from a file that contains edges as space separated pairs of node indexes in the format: `source_index  target_index`. The labels are loaded from a file that contains the label of each node as `node_index node_label` in each line.  

__OUTPUT FILE FORMAT__

Output is file containing the summary graph.

__COMPILATION__ (LINUX)

Dependencies: `boost` library must be installed

Command line:  `make`

__EXECUTION__
		      	 
Command line: `./sum [OPTIONS]`

__ARGUMENTS & OPTIONS__

ARGUMENT | TYPE | DEFAULT VALUE
-------- | ------ | -------
`--graph_file` | path to file containing edge-list (see above) | `"../../../../graphs/HomoSapiens/adj.txt"`
`--outfile` | (embeddings) | `"../../embeddings/HomoSapiens_embed.txt"`
`--dimension` | Dimension of embedding space (vector length) | `100`
`--walk_length` | Maximum length of walks (similarity powers) considered | `10`
`--lambda` | l-2 regularization parameter | `1.0e-3`
`--splits` | Number of splits | `2`
`--edges_per_split` | Number of positive edges removed at every split | `1000`     
`--fit_with` | Select model according to which proximity parameters are learnt: 1) Logistic regression (`logistic`), 2) Least squares (`ls`), 3) SVMs `svm`, 4) Chose single best proximity/walk-length (`single best`)  | `single_best`

