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
#include <unordered_map>
#include <cassert>
#include <cmath>

#include "node.h"
#include "utils.h"

using namespace std;



/*------------------------------------------------------------------------

            NODE METHODS

-------------------------------------------------------------------------*/

/*
Use the inner structure (+ maybe the connections) to decide on the supenode glyph
*/
void node::decide_glyph(void)
{
    if(this->simple_node)
    {
        this->glyph = NONE;
        return;
    } 

    /*
    Get pair() of in- and out- degrees for every node of inner graph (excluding self-loops)
    */
    unordered_map<int,pair<size_t,size_t>> deg;
    for(auto it = this->nodeset.cbegin(); it != this->nodeset.cend(); ++it)
        deg.insert(pair<int,pair<size_t,size_t>>( *it, pair<size_t,size_t>(0,0)));

    size_t total_edges = 0 ;
    for(auto it = this->ingraph.cbegin(); it != this->ingraph.cend(); ++it)
    {
        for(auto it2 = it->second.cbegin(); it2 != it->second.cend(); ++it2)
        {
            if(it->first != it2->first)
            {
                deg[it->first].second ++;
                deg[it2->first].first ++;
                total_edges ++;
            }
        }
                        
    }


    // cout << "\n\n\nnode : ";
    // for(auto it = deg.cbegin(); it != deg.cend(); ++it)
    //     cout << " " << it->first;
    // cout << "\nin-deg : ";
    // for(auto it = deg.cbegin(); it != deg.cend(); ++it)
    //     cout << " " << it->second.first;
    // cout << "\nout-deg : ";
    // for(auto it = deg.cbegin(); it != deg.cend(); ++it)
    //     cout << " " << it->second.second;


    /* 
        Check if it is a near-CLIQUE
    */

    if(total_edges >= this->size *(this->size-1)/2)
    {
        this->glyph = CLIQUE;
        return;
    }

    /*
        Check for INSTAR or OUTSTAR
    */
    size_t max_out = 0, max_in = 0;
    int in_hub = 0, out_hub = 0;
    for(auto it = deg.cbegin(); it != deg.cend(); ++it)
    {
        if(it->second.first > max_in)
        {
            max_in = it->second.first;
            in_hub = it->first;
        }
        if(it->second.second > max_out)
        {
            max_out = it->second.second;
            out_hub = it->first;
        }
    }

    double instar_cost_approx = (  abs((double)(this->size - 1 - max_in)) +  abs((double)(total_edges - max_in)) ) / (double) (total_edges);
    double outstar_cost_approx = (  abs((double)(this->size - 1 - max_out)) +  abs((double)(total_edges - max_out)) ) / (double) (total_edges);

    // cout << "\n\ninstar cost: " << instar_cost_approx ;
    // cout << "\n\noutstar cost: " << outstar_cost_approx ;    

    double thres = 1.0;
    if(instar_cost_approx < thres || outstar_cost_approx < thres )
    {
        if(instar_cost_approx < outstar_cost_approx )
        {
            this->glyph = INSTAR;
            this->hub = in_hub;
        }
        else
        {
            this->glyph = OUTSTAR;
            this->hub = out_hub;
        }
        return;
    }

    /*
        If not a star then DISConnected
    */
    this->glyph = DISC;
    return;
}



/*
Obtain representative glyph and self-loop multiplicities  
using dichotomous search on the [min,max] range
*/
void node::get_rep_mult(void)
{
    /*
        Extract self-loop and glyph multiplicities
    */
    vector<size_t> sl_mults;
    vector<size_t> gl_mults;
    for(auto it = this->ingraph.cbegin(); it != this->ingraph.cend(); ++it)
    {
        for(auto it2 = it->second.cbegin(); it2 != it->second.cend(); ++it2)
        {
            if(it->first != it2->first)
                gl_mults.push_back(it2->second);
            else
                sl_mults.push_back(it2->second);                            
        }
                        
    }   

    /*
        Call dichotomous on self-loop multiplicities 
    */
    if(sl_mults.size())
    {
        pair<size_t, double> sol1 = dichotomous_search(sl_mults); 

        this->sl_mult = sol1.first;
        this->sl_mult_cost = sol1.second;
    }

    /*
        Call dichotomous on glyph multiplicities 
    */
    if(gl_mults.size())
    {
        pair<size_t, double> sol2 = dichotomous_search(gl_mults); 

        this->nd_mult = sol2.first;
        this->nd_mult_cost = sol2.second;
    }   

}







/*
Obtain or update the encoding cost for rewiring edges of this glyph ( self loop seperate to glyphs)  
*/
void node::get_rewiring_cost()
{

    if(this->glyph == NONE)
    {
        assert(this->size == 1);
        return ;
    }

    /*
    Get pair() of in- and out- degrees for every node of inner graph (excluding self-loops)
    Also measure total true edges and total self loops
    */
    unordered_map<int,pair<size_t,size_t>> deg;
    for(auto it = this->nodeset.cbegin(); it != this->nodeset.cend(); ++it)
        deg.insert(pair<int,pair<size_t,size_t>>( *it, pair<size_t,size_t>(0,0)));

    size_t total_edges = 0;
    size_t total_self_loops = 0;
    for(auto it = this->ingraph.cbegin(); it != this->ingraph.cend(); ++it)
    {
        for(auto it2 = it->second.cbegin(); it2 != it->second.cend(); ++it2)
        {
            if(it->first != it2->first)
            {
                deg[it->first].second ++;
                deg[it2->first].first ++;
                total_edges ++;
            }
            else
            {
                total_self_loops ++;
            }
            
        }
                        
    }
   
    /*
    First get self-loop related costs
    */
    this->sl_cost = binom_cost( this->size, this->size - total_self_loops);
    this->non_sl_cost = binom_cost( this->size, total_self_loops );
 

    /*
    Continue wih cost per glyph
    */
    if(this->glyph == CLIQUE)
    {
        this->glyph_cost = binom_cost( this->size * (this->size - 1) , this->size * (this->size - 1) - total_edges);
        return ;
    }

    if(this->glyph == DISC)
    {
        this->glyph_cost = binom_cost( this->size * (this->size - 1) , total_edges);
        return ;
    }

    if(this->glyph == INSTAR)
    {
        size_t indeg = deg[this->hub].first;
        this->glyph_cost = 0.0;
        /*
            Positive corrections
        */
        this->glyph_cost += binom_cost( this->size * (this->size - 1) - (this->size - 1), total_edges - indeg) 
                            - (total_edges - indeg > 0 ? log2((double)(total_edges - indeg)) : 0) - 1;
        //cout << "\n\npos_cor " << this->glyph_cost << "\n";

        /*
            Negative corrections
        */
        this->glyph_cost += binom_cost( this->size - 1, this->size - 1 - indeg) 
                            - (this->size - 1 - indeg > 0 ? log2((double)(this->size - 1 - indeg)) : 0) - 1;
        // cout << "\n\nneg_cor " << binom_cost( this->size - 1, this->size - 1 - indeg) << "\n\n";

       // cout << "\n\nnode_size " << this->size << ", total_edges " << total_edges << ", indeg " << indeg << "\n";
        return ;
    }

    if(this->glyph == OUTSTAR)
    {
        size_t outdeg = deg[this->hub].second;
        this->glyph_cost = 0.0;
        /*
            Positive corrections
        */
        this->glyph_cost += binom_cost( this->size * (this->size - 1) - (this->size - 1), total_edges - outdeg)
                        - (total_edges - outdeg > 0 ? log2((double)(total_edges - outdeg)) : 0) - 1;
        /*
            Negative corrections
        */
        this->glyph_cost += binom_cost( this->size - 1, this->size - 1 - outdeg)
                         - (this->size - 1 - outdeg > 0 ? log2((double)(this->size - 1 - outdeg)) : 0) - 1;
        return ;
    }

}


/*
  (Re)-Compute and sort edge_delta vector
*/
void node::update_edge_deltas(void)
{
    /* Remove previous values */
    this->sedge_deltas.clear();

    /* Place new ones */
    this->sedge_deltas.reserve(this->sedges.size()+1);

    double sl_delta = this->sl_cost - this->non_sl_cost;
    sl_delta += this->sl_mult_cost - this->non_sl_mult_cost;        
    //if(sl_delta < 0.0)
        this->sedge_deltas.push_back( pair<int, double>(SELF_ID, sl_delta));


    for (auto it= this->sedges.cbegin(); it != this->sedges.cend(); ++it)
    {
        double delta = it->second.rewiring_cost - it->second.non_rewiring_cost;
        delta += it->second.rep_mult_cost - it->second.non_rep_mult_cost;        
       // if(delta < 0.0)
            this->sedge_deltas.push_back( pair<int, double>( it->first, delta));        
    }
    

    /* Sort in ncreasing order wrt deltas*/
    sort(this->sedge_deltas.begin(), this->sedge_deltas.end(), 
         [](const pair<int, double>& lhs, const pair<int, double>& rhs) {return lhs.second < rhs.second; } );
}

/*
    (Re)-evaluate edge activation
    NOTE: edge_deltas MUST be sorted in increasing order
*/
void node::activate_edges( size_t nnodes)
{
    /* First deactivate all super-edges */
    this->self_loop = false;
    for(auto it= this->sedges.begin(); it != this->sedges.end(); ++it)
        it->second.active = false;

    /* Procceed with activation */
    if(this->sedge_deltas.empty())
        return;

#if(false)
    /*
    TODO: Carefully look at edge cases..
    */
    // if(this->sedge_deltas.size() == 1)
    // {
        // if(this->sedge_deltas[0].second + log2((double) nnodes) < 1.0)
        // {
        //     this->sedges.
        // }    

    // }

    double min_cost = 1.0;
    size_t min_at = 0;
    double sum_at = 1.0;

    for(size_t i=0; i< min(this->sedge_deltas.size(), nnodes/2) ; i++)
    {
        /* Update cost */
        sum_at += this->sedge_deltas[i].second;
        sum_at += log2( (double) ( nnodes - (i+1) )) - log2( (double) (i+2));

        /* Track minimum */
        if(sum_at < min_cost)
        {
            min_cost = sum_at;
            min_at = i+1;
        }
    }

    /* Activate first min_at super-edges */
    for(size_t i=0; i< min_at; i++)
    {
        int sedge_to = this->sedge_deltas[i].first;
        if(sedge_to == SELF_ID)
            this->self_loop = true;
        else
            this->sedges[sedge_to].active = true;
    }
#else

    /* Activate all negative deltas */
    for(size_t i=0; i< this->sedge_deltas.size(); i++)
    {
        if(this->sedge_deltas[i].second < 0.0)
        {
            int sedge_to = this->sedge_deltas[i].first;
            if(sedge_to == SELF_ID)
                this->self_loop = true;
            else
                this->sedges[sedge_to].active = true;
        }
    }

#endif

}





















