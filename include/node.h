#ifndef NODE_H
#define NODE_H

#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include <unordered_set>
#include <unordered_map>
#include <set>



#include "edge.h"
#include "utils.h"

/*------------------------------------------------------------------------

                                    NODE

-------------------------------------------------------------------------*/



// Super-node data structure
class node
{
  public:
    /*
       General Supernode Summary Information
    */
    bool simple_node;
    int type;
    size_t size;
    size_t nd_mult;
    bool self_loop;
    size_t sl_mult;
    enum glyph_type glyph;

    /* Costs */
    double nd_mult_cost;
    double sl_mult_cost;
    double non_sl_mult_cost;
    double glyph_cost;
    double sl_cost;
    double non_sl_cost;

    /* INNER-STRUCTURE

       1) hub: denotes which (original) node ID that plays the role of hub 

       2) ingraph: key is original node id contained in super-node, value is list of nodes it points at and multiplicities

    */ 
    int hub; 
    std::set<int> nodeset;
    std::unordered_map<int,std::unordered_map<int, size_t>> ingraph;

    /*
    Outgoing Superedges (key: destination (super) node, value: super-edge object)
    */
    std::vector<std::pair<int, double>> sedge_deltas;
    std::unordered_map<int, edge> sedges;

    /*
     IDs of super-nodes it is pointed by
    */
    std::unordered_set<int> pointed_by;

  public:
    /*
    Default constructor as  simple node
    */
    node() : simple_node{true}, type{0}, size{1}, nd_mult{1}, self_loop{false},
             sl_mult{1}, glyph{NONE}, nd_mult_cost{0.0}, sl_mult_cost{1.0},
             non_sl_mult_cost{1.0}, glyph_cost{0.0}, sl_cost{1.0}, non_sl_cost{0.0}, hub{0}
    {
        ingraph.reserve(1);
        sedge_deltas.reserve(1);
        sedges.reserve(1);
        pointed_by.reserve(1);
    }

    /*
    Use the inner structure (+ maybe the connections) to decide on the supenode glyph
    */
    void decide_glyph();  

    /*
    Obtain representative node and self-loop multiplicities multiplicities 
    using dichotomous search on the [min,max] range
    */
    void get_rep_mult();

    /*
    Obtain or update the encoding cost for rewiring edges of this glyph ( self loop seperate to glyphs)  
    */
    void get_rewiring_cost();

    /*
      (Re)-Compute and sort edge_delta vector
    */
    void update_edge_deltas();

    /*
        (Re)-evaluate edge activation
        NOTE: edge_deltas MUST be sorted in increasing order
    */
    void activate_edges( size_t nnodes);

};


// #include "utils.h"

#endif