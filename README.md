__INPUT FILES FORMAT__

__EMB_min__ loads the graph in __edge list__ format from a file that contains edges as space separated pairs of node indexes in the format: `source_index  target_index`. The labels are loaded from a file that contains the label of each node as `node_index node_label` in each line.  

__OUTPUT FILE FORMAT__

Output is file containing the summary graph. First, the description of each super-node in every line, followed
by the description of each super-edge in every line. See bellow for line formats.


SUPER-NODES 
`ID  Size  Type  Glyph Glyph-rep-mult Self-loop SL-rep-mult`

SUPER-EDGES
`Source-ID  Dest-ID Rep-mult`


__COMPILATION__ (LINUX)

Dependencies: `boost` library must be installed

Command line:  `make`

__EXECUTION__
		      	 
Command line: `./sum [OPTIONS]`

__ARGUMENTS & OPTIONS__

ARGUMENT | DESCRIPTION
-------- | ------ | -------
`--graph` | Path to file containing edge-list (see above) 
`--save-to` | Path to summary graph output
`--bands` | Set number (>=1) of lsh bands
`--rows` | Set number (>=1) of lsh rows
`--K` | Maximum number of candidate sets to be considered
`--niters` | Number of experiments
