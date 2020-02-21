#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <queue>
#include <numeric>

#include "edge.h"
#include "node.h"
// #include "utils.h"

/*------------------------------------------------------------------------

                                    GRAPH 

-------------------------------------------------------------------------*/

// Graph with adjacency list representation for effiecient random walks
// Can also be used for APPR with push() prodedure for fast personalized PageRank
class graph
{
  public:
    const static int init_size = 100;

    // Adjacency list representation
    std::set<int> ids;
    std::unordered_map<int,node> snode;
    size_t nnodes_full;
    int nnodes;
    size_t ntypes;
    bool no_mults;

  public:
    //Default constr.
    graph() : nnodes_full{0}, nnodes{0}, ntypes{0}, no_mults{false}
    {
        snode.reserve(init_size);
    }

    //Constr. from file
    graph(std::string graph_file, std::string types_file, bool no_mults, bool no_types);

    /*
     Observers
     */
    void display();
    /*  Obtain  encoding cost of summary  */
    std::pair<double,double> get_cost();

    /*  Obtain  simple copy of forward super-graph topology  */
    std::unordered_map<int,std::vector<int>> get_fwd_graph();

    /*  Obtain  simple copy of inverse super-graph topology  */
    std::unordered_map<int,std::vector<int>> get_rev_graph();


    /*
     Modifiers
    */

    /* Merge set of super-nodes/nodes (return new super-node id) */ 
    int merge_snodes(std::vector<int>& mrg_ids);

    /*  Expand supernodes to the simple nodes and edges they contain */
    void expand_snodes(std::vector<int>& exp_ids);

    /* Cleanup empty supernodes */
    void cleanup();

    /*
      Storage
    */

    /* Save summary graph to file*/
    void save_as(std::string path);

};




#endif
