#ifndef EDGE_H
#define EDGE_H

#include <string>
#include <vector>
#include <memory>
#include <cstdint>


#include "utils.h"

/*------------------------------------------------------------------------

                                    EDGE

-------------------------------------------------------------------------*/

// Super-edge data structure
class edge
{
  public:
    /*
       Super-edge Information
    */
    bool simple_edge;
    size_t size;
    size_t rep_mult;
    double rep_mult_cost;
    double non_rep_mult_cost;
    double rewiring_cost;
    double non_rewiring_cost;

   /* 
      Whether it is active in summary (ignore for simple-edge)
   */

    bool active;   

  public:
    /* 
      CONTENTS : List of simple edges as the triplet (source, dest, mult)
    */ 

    std::vector<int> src_node;
    std::vector<int> dest_node;
    std::vector<size_t> mults;

    /*
    Default constructor as  simple node
    */
    edge() : simple_edge{true}, size{1}, rep_mult{1}, rep_mult_cost{0.0},
             non_rep_mult_cost{0.0}, rewiring_cost{0.0}, non_rewiring_cost{1.0}, active{true}
    {
      src_node.reserve(2);
      dest_node.reserve(2);
      mults.reserve(2);
    }

    /*
    Obtain representative superdge multiplicity 
    using dichotomous search on the [min,max] range
    */
    void get_rep_mult();

    /*
    Obtain or update the encoding cost for rewiring edges of this superedge (if active)  
    */
    void get_rewiring_cost(enum glyph_type glyph_from, enum glyph_type glyph_to, size_t size_from,
                           size_t size_to, int hub_from, int hub_to, std::set<int>& nodeset_from, std::set<int>& nodeset_to  );

};

#endif