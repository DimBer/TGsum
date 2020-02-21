#include <cstdlib>
#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <random>
#include <iomanip>
#include <unordered_set>
#include <cmath>
#include <functional>

#include "node.h"
#include "edge.h"
#include "utils.h"


using namespace std;

/*------------------------------------------------------------------------

            EDGE METHODS

-------------------------------------------------------------------------*/

/*
Obtain representative super-edge  multiplicity 
using dichotomous search on the [min,max] range
*/
void edge::get_rep_mult()
{
    if(this->simple_edge || this->size == 1 )
    {
       if(this->mults.size()) this->rep_mult = this->mults[0]; else this->rep_mult = 1; 
       this->rep_mult_cost = 2.0 * log2( (double) this->rep_mult) + 1.0;
       this->non_rep_mult_cost = this->rep_mult_cost;
       return;
    }    

    /*
        Extract edge multiplicities
    */
    this->non_rep_mult_cost = 0.0;
    for(auto& m : this->mults)
        this->non_rep_mult_cost += 2.0 * log2( (double) m) + 1.0;

    /*
        Call dichotomous on multiplicities 
    */
    pair<size_t, double> sol = dichotomous_search(this->mults); 

    this->rep_mult = sol.first;
    this->rep_mult_cost = sol.second;
}


/*
Obtain or update the encoding cost for rewiring edges of this superedge (if active)  
REQUIRES: Info on size and glyph of adjacent supernodes
*/
void edge::get_rewiring_cost(enum glyph_type glyph_from, enum glyph_type glyph_to, size_t size_from,
                             size_t size_to, int hub_from, int hub_to, set<int>& nodeset_from_, set<int>& nodeset_to_ )
{

    unordered_set<int> nodeset_from(nodeset_from_.cbegin(),nodeset_from_.cend());
    unordered_set<int> nodeset_to(nodeset_to_.cbegin(),nodeset_to_.cend());

    if(glyph_from == NONE && glyph_to == NONE )
        return ;

    this->non_rewiring_cost = binom_cost( size_from*size_to, this->size) - log2((double) this->size);


    size_t pos_cor = 0, neg_cor = 0, expand_size = 0; 


    if(glyph_from == CLIQUE || glyph_from == DISC || glyph_from == NONE)
    {
        if(glyph_to == CLIQUE || glyph_to == DISC || glyph_to == NONE)
        {
            /*
             Expands  all-to-all
            */
            pos_cor = 0;
            neg_cor = size_from*size_to - this->size;
            expand_size = size_from*size_to;
        }
        else 
        {
            /*
             Expands  all-to-hub
            */
            pos_cor = 0;
            for(size_t i=0; i<this->size; i++)
            { 
                if(this->dest_node[i] != hub_to && this->dest_node[i] != SELF_ID )
                    pos_cor++;    
                nodeset_from.erase(this->src_node[i]);
            }   
            neg_cor = nodeset_from.size();
            expand_size = size_from;        
        }    
    }    
    else 
    {
        if(glyph_to == CLIQUE || glyph_to == DISC || glyph_to == NONE)
        {
            /*
             Expands  hub-to-all
            */
            pos_cor = 0;
            for(size_t i=0; i<this->size; i++)
            { 
                if(this->src_node[i] != hub_from && this->src_node[i] != SELF_ID )
                    pos_cor++;    
                nodeset_to.erase(this->dest_node[i]);
            }   
            neg_cor = nodeset_to.size();
            expand_size = size_to;      
            // cout << "hub from " << hub_from << ", pos_cor " << pos_cor << ", neg_cor " << neg_cor <<  ", edge size " << this->size  << "\n";
        }
        else 
        {
            /*
             Expands  hub-to-hub
            */
            pos_cor = 0;
            bool h2h_found = false;
            for(size_t i=0; i< this->size; i++)
            {
                if( (this->src_node[i] == hub_from) && (this->dest_node[i] == hub_to))
                    h2h_found = true;
                else
                    pos_cor++;                
            }
            neg_cor = h2h_found ? 0 : 1 ;
            expand_size = 1;                
        }       
    }

    
    this->rewiring_cost = 0.0;
    this->rewiring_cost += binom_cost( size_from*size_to - expand_size , pos_cor);
    this->rewiring_cost += binom_cost( expand_size , neg_cor);
    // cout << "\n\npos cor " <<  pos_cor << " , neg cor " << neg_cor <<", expands " << expand_size << "\n";
    // cout << "\n\rewire cost " <<  rewiring_cost <<  "\n\n";
}