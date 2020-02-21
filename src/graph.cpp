#include <cstdlib>
#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <random>
#include <cassert>
#include <set>
#include <sys/time.h>

#include "edge.h"
#include "node.h"
#include "utils.h"
#include "graph.h"

#define EXPORT_SELF_LOOPS true

#define PROFILE_MERGE false

using namespace std;

/*
Predicate for using max_element
*/
static bool pred(const std::pair<int,node>& lhs, const std::pair<int,node>& rhs)
{
    return lhs.first < rhs.first;
}


/*------------------------------------------------------------------------

                 GRAPH METHODS

-------------------------------------------------------------------------*/

/*
Constructor of graph object
Loads edgelist file and builds adjacency list
*/
graph::graph(string graph_file, string types_file, bool no_mults, bool no_types)
{
    this->no_mults = no_mults;

    this->nnodes = 0;
    
    /*
    Load graph structure
    */
    cout << "Loading graph from: " << graph_file << "\n";
    ifstream myfile(graph_file);
    int node_a{0}, node_b{0};
    if (myfile.is_open())
    {
        while (!myfile.eof())
        {
            node_a = 0;
            node_b = 0;
            myfile >> node_a;
            myfile >> node_b;

            if(node_a == 0 && node_b == 0 )
                continue;

            /*
             Add new node(s) if it does not exist
            */
            this->ids.insert(node_a);
            this->ids.insert(node_b);

            /*
             Add  node_a if it doesnt exist yet
            */
            if(!this->snode.count(node_a))
            {
                this->snode.insert(pair<int,node>(node_a, node()));
                this->nnodes ++;
            }
                

            if(node_a != node_b)
            {
                /*
                Record edge as simple superedge
                */
                if(!this->snode[node_a].sedges.count(node_b))
                {
                    this->snode[node_a].sedges.insert(pair<int, edge>(node_b, edge()));
                    this->snode[node_a].sedges[node_b].src_node.push_back(node_a);
                    this->snode[node_a].sedges[node_b].dest_node.push_back(node_b);
                    this->snode[node_a].sedges[node_b].mults.push_back(1);

                    /* Initialize super-edge delta */
                    this->snode[node_a].sedge_deltas.reserve(this->snode[node_a].sedge_deltas.size() + 1);
                    this->snode[node_a].sedge_deltas.push_back(pair<int, double>(node_b, -1.0) );                    
                }    
                else
                {
                    for(size_t i=0; i<this->snode[node_a].sedges[node_b].size; i++)
                        if(this->snode[node_a].sedges[node_b].src_node[i] == node_a && this->snode[node_a].sedges[node_b].dest_node[i] == node_b )
                        {
                            /* Increase multiplicities and related initial costs */
                            this->snode[node_a].sedges[node_b].mults[i] += 1;
                            this->snode[node_a].sedges[node_b].rep_mult ++;
                            this->snode[node_a].sedges[node_b].rep_mult_cost = 2.0*log2( (double) this->snode[node_a].sedges[node_b].rep_mult)  +1.0;
                            this->snode[node_a].sedges[node_b].non_rep_mult_cost = this->snode[node_a].sedges[node_b].rep_mult_cost;
                            // cout << "\n here from" << node_a << " to " << node_b << " cost " << this->snode[node_a].sedges[node_b].rep_mult_cost;
                        }
                            
                }

                /*
                 Add node_b to the pointed_by of node_a
                */
                if(!this->snode.count(node_b))
                {
                    this->snode.insert(pair<int,node>(node_b, node()));
                    this->nnodes ++;
                }
                this->snode[node_b].pointed_by.insert(node_a);

            }
            else
            {
                /*
                Self loops are recorded as part of inner structure
                */
                this->snode[node_a].self_loop = true;
                if(!this->snode[node_a].ingraph.count(node_a))
                {
                    this->snode[node_a].ingraph.insert(pair<int,unordered_map<int,size_t>>(node_a, unordered_map<int,size_t>()));
                    this->snode[node_a].ingraph[node_a].insert(pair<int,size_t>(node_a,1));

                    /* Initialize self-loop  delta */
                    this->snode[node_a].self_loop = true;
                    this->snode[node_a].sl_cost = 0.0;
                    this->snode[node_a].non_sl_cost = 1.0;
                    this->snode[node_a].sedge_deltas.reserve(this->snode[node_a].sedge_deltas.size() + 1);
                    this->snode[node_a].sedge_deltas.push_back(pair<int, double>(SELF_ID, -1.0) );
                }    
                else
                {
                    if(!this->snode[node_a].ingraph[node_a].count(node_a))
                    {
                        this->snode[node_a].ingraph[node_a].insert(pair<int,size_t>(node_a,1));
                    }
                    else
                    {
                        this->snode[node_a].ingraph[node_a].at(node_a) += 1;
                        this->snode[node_a].sl_mult++;
                        this->snode[node_a].sl_mult_cost = 2.0*log2( (double) this->snode[node_a].sl_mult) + 1.0;
                        this->snode[node_a].non_sl_mult_cost = this->snode[node_a].sl_mult_cost;
                    }
                }                
            }
            

            /*
             Initialize node sets of super-node a and b wit them-selves
            */
            this->snode[node_a].nodeset.insert(node_a);
            this->snode[node_b].nodeset.insert(node_b);
        
        }
        myfile.close();
    }
    else
        cout << "Unable to open file";

    /* Get rep mults and correct deltas for all super-nodes */
    for(auto nodeit = this->snode.begin(); nodeit != this->snode.end(); ++nodeit)
    {
        for(auto edgeit = nodeit->second.sedges.begin(); edgeit != nodeit->second.sedges.end(); ++edgeit)
            edgeit->second.get_rep_mult();
        nodeit->second.update_edge_deltas();
        nodeit->second.activate_edges(this->nnodes);
    }

    /*
    Load types
    */
    set<int> type_set;
    if(!no_types)
    {
        cout << "Loading types from: " << types_file << "\n";
        ifstream myfile2(types_file);
        int node_id{0}, type{0};
        if (myfile2.is_open())
        {
            while (!myfile2.eof())
            {
                node_id = 0;
                type = 0;
                myfile2 >> node_id;
                myfile2 >> type;

                if(node_id == 0 && type == 0 )
                    continue;

                if(!this->snode.count(node_id))
                {
                    cout << "WARNING: Node %d appears in type file but not in graph file\n";        
                }
                else
                {
                    this->snode[node_id].type = type;
                    type_set.insert(type);
                }
                

            }
            myfile2.close();
        }
        else
            cout << "Unable to open file";
    }
    else
    {
        type_set.insert(1);
        for(auto it= this->snode.begin(); it != this->snode.end(); ++it)
            it->second.type = 1;
    }

    this->ntypes = type_set.size();
    this->nnodes_full = this->nnodes;
}





/*
Display the graph
*/
void graph::display(void)
{
    cout << "DISPLAYING GRAPH\n\n";
    cout << "# supernodes: " << this->snode.size() << "\n\n";

    for(auto nodeit = this->snode.cbegin(); nodeit != this->snode.cend(); ++nodeit )
    {
        cout << "\nSuper-node ID: " << nodeit->first ; 
        cout << "\nType: " << nodeit->second.type; 
        cout << "\nSize: " << nodeit->second.size; 
        cout << "\nGlyph multiplicity: " << nodeit->second.nd_mult;     
        cout << "\nSelf loop: ";
        if(nodeit->second.self_loop) cout << "YES"; else cout<< "NO";         
        cout << "\nSelf loop mult.: " << nodeit->second.sl_mult;
        cout << "\nGlyph: ";
        if(nodeit->second.glyph == DISC )
            cout << "DISCONNECTED";
        else if(nodeit->second.glyph == CLIQUE )
            cout << "CLIQUE";            
        else if(nodeit->second.glyph == INSTAR )
            cout << "INSTAR with hub: " << nodeit->second.hub;            
        else if(nodeit->second.glyph == OUTSTAR )
            cout << "OUTSTAR with hub: " << nodeit->second.hub;            
        else
            cout << "NONE";        

        cout << "\n\nNodes Contained: ";
        for(auto innodeit = nodeit->second.nodeset.cbegin(); innodeit != nodeit->second.nodeset.cend(); ++innodeit)
            cout << " " << *innodeit ;

        cout << "\n\nPoints to: ";
        for(auto edgeit = nodeit->second.sedges.cbegin(); edgeit != nodeit->second.sedges.cend(); ++edgeit)
            cout << " " << edgeit->first << "(" << edgeit->second.rep_mult << "/" << edgeit->second.rewiring_cost - edgeit->second.non_rewiring_cost 
                 << "/" << (int)edgeit->second.active << ")";

        cout << "\n\nPointed by: ";
        for(auto inedgeit = nodeit->second.pointed_by.cbegin(); inedgeit != nodeit->second.pointed_by.cend(); ++inedgeit)
            cout << " " << *inedgeit ;

        cout<< "\n\n";

    }

    pair<double,double> costs = this->get_cost();
    cout << "COST = [ summ: " << costs.first << ", corr: " << costs.second << ", total:  " << costs.first + costs.second << "]\n\n";

}



/************************************************************************************************
    Merging subset of super-nodes
*************************************************************************************************/
int graph::merge_snodes(vector<int>& mrg_ids)
{
    if(mrg_ids.size() <= 1 )
        return 0;

#if PROFILE_MERGE
    /* Runtime profiling */
    struct timeval start = {.tv_sec = 0, .tv_usec = 0 };
    struct timeval stop = {.tv_sec = 0, .tv_usec = 0 };
#endif

    /**********************************************************
    PHASE 1
    ***********************************************************/

    /*
        Make sure all super-nodes have the same type
    */
        for(auto& id : mrg_ids)
            assert( this->snode[id].type == this->snode[mrg_ids[0]].type );

    /*
        Create the new supernode
    */
        int new_node = max_element(this->snode.cbegin(),this->snode.cend(),pred)->first + 1;
        this->snode.insert(pair<int,node>(new_node, node()));
        this->ids.insert(new_node);
        this->nnodes++;

    /**********************************************************
    PHASE 2
    ***********************************************************/

#if PROFILE_MERGE   
    gettimeofday( &start, NULL);
#endif

    /*
        Populate the new supernode
    */
    this->snode[new_node].simple_node = false;
    this->snode[new_node].type = this->snode[mrg_ids[0]].type;
    this->snode[new_node].size = 0;
    for(auto& id : mrg_ids)
    {
        /*
            Merging nodesets and addings sizes
        */
        this->snode[new_node].nodeset.insert(this->snode[id].nodeset.begin(),this->snode[id].nodeset.end());
        this->snode[new_node].size += this->snode[id].size; 

        /*
            Merge inner structures (intersection should be zero)
        */
        this->snode[new_node].ingraph.insert(this->snode[id].ingraph.begin(),this->snode[id].ingraph.end());

        /*
            Merging the pointed_by lists (union of sets)
        */
        this->snode[new_node].pointed_by.insert(this->snode[id].pointed_by.begin(),this->snode[id].pointed_by.end());
    }

    /*
        Remove merged nodes from merged pointed_by list
    */
    for(auto& id : mrg_ids)
        this->snode[new_node].pointed_by.erase(id);

    /*
        Merge outgoing superedges
    */
    set<int> mrg_ids_set;
    mrg_ids_set.insert(mrg_ids.cbegin(),mrg_ids.cend());
    for(auto& id : mrg_ids)
    {
        /*
            Go through every super-edge
        */
        for(auto sedgeit = this->snode[id].sedges.cbegin(); sedgeit != this->snode[id].sedges.cend(); ++sedgeit)
        {     
            if(mrg_ids_set.count(sedgeit->first)>0)
            {
               // cout << "Redirecting to ingraph\n";
                /*
                    In-between superedges are redirected to the new ingraph
                */
                for(size_t i=0; i<sedgeit->second.size; i++)
                {
                    if(!this->snode[new_node].ingraph.count(sedgeit->second.src_node[i]))
                        this->snode[new_node].ingraph.insert(pair<int,unordered_map<int,size_t>>(sedgeit->second.src_node[i], unordered_map<int,size_t>()));
                    this->snode[new_node].ingraph[sedgeit->second.src_node[i]].insert( pair<int, size_t>(sedgeit->second.dest_node[i], sedgeit->second.mults[i]));
                }
            }
            else
            {
                /*
                    Merge the TRULY outgoing super-edges 
                */
                if(!this->snode[new_node].sedges.count(sedgeit->first))
                {
                    this->snode[new_node].sedges.insert( pair<int, edge>(sedgeit->first, sedgeit->second ));
                }
                else
                {
                    //cout << "Merging super-edge\n";

                    this->snode[new_node].sedges[sedgeit->first].size += sedgeit->second.size;

                    auto& new_src =  this->snode[new_node].sedges[sedgeit->first].src_node;
                    auto& new_dest =  this->snode[new_node].sedges[sedgeit->first].dest_node;
                    auto& new_mults =  this->snode[new_node].sedges[sedgeit->first].mults;                

                    new_src.reserve(new_src.size() + sedgeit->second.src_node.size());
                    new_dest.reserve(new_dest.size() + sedgeit->second.dest_node.size());
                    new_mults.reserve(new_mults.size() + sedgeit->second.mults.size());

                    for(size_t i=0; i< sedgeit->second.size; i++ )
                    {
                        new_src.push_back(sedgeit->second.src_node[i]);
                        new_dest.push_back(sedgeit->second.dest_node[i]);
                        new_mults.push_back(sedgeit->second.mults[i]);
                    }
                }
            }
            
        }
    }

#if PROFILE_MERGE
    gettimeofday( &stop, NULL);
    auto dur =  stop.tv_usec - start.tv_usec;
    cout << "\n\nphase 1: " << (size_t) dur;
#endif

    /*
    NOTE: At this point the merged super-nodes are ``coexisting'' 
        with their merged version 
    */


    /**********************************************************
    PHASE 3
    ***********************************************************/
#if PROFILE_MERGE   
    gettimeofday( &start, NULL);
#endif

    /*
    Find one hop neighborhoods of merged nodes (excluding themselves)
    */
    set<int> neigh_set(this->snode[new_node].pointed_by.begin(), this->snode[new_node].pointed_by.end());
    for(auto edgeit = this->snode[new_node].sedges.cbegin(); edgeit != this->snode[new_node].sedges.cend(); ++edgeit)
        neigh_set.insert(edgeit->first);

    /*
    Assert that neighbors and merged nodes have zero intersection
    */
    vector<int> neigh(neigh_set.cbegin(),neigh_set.cend()); 
    sort(mrg_ids.begin(),mrg_ids.end(), greater<int>());
    sort(neigh.begin(),neigh.end(), greater<int>());
    vector<int> intrs(max(mrg_ids.size(), neigh.size()));
    auto it = set_intersection(mrg_ids.cbegin(),mrg_ids.cend(), neigh.cbegin(), neigh.cend(), intrs.begin());
    intrs.resize(it-intrs.begin());

    assert(intrs.size() == 0 );

    /*
    Continue to updates..
    */

    for(auto& this_neigh : neigh)
    {
        bool is_pointed = false;
        for(auto& this_mrg : mrg_ids)
        {
            /*
                Update neighbor super-node FROM references
            */
            if(this->snode[this_neigh].pointed_by.erase(this_mrg) > 0) is_pointed = true;

            /*
                Update neighbor super-node TO references
            */

            if(this->snode[this_neigh].sedges.count(this_mrg) > 0)
            {
                if(this->snode[this_neigh].sedges.count(new_node) == 0)
                {
                    this->snode[this_neigh].sedges.insert( pair<int, edge>( new_node, this->snode[this_neigh].sedges[this_mrg]));                
                }
                else
                {
                    this->snode[this_neigh].sedges[new_node].size += this->snode[this_neigh].sedges[this_mrg].size;

                    auto& new_src =  this->snode[this_neigh].sedges[new_node].src_node ;
                    auto& new_dest =  this->snode[this_neigh].sedges[new_node].dest_node;
                    auto& new_mults =  this->snode[this_neigh].sedges[new_node].mults;                

                    new_src.reserve(new_src.size() + this->snode[this_neigh].sedges[this_mrg].src_node.size());
                    new_dest.reserve(new_dest.size() + this->snode[this_neigh].sedges[this_mrg].dest_node.size());
                    new_mults.reserve(new_mults.size() + this->snode[this_neigh].sedges[this_mrg].mults.size());

                    for(size_t i=0; i< this->snode[this_neigh].sedges[this_mrg].size; i++ )
                    {
                        new_src.push_back(this->snode[this_neigh].sedges[this_mrg].src_node[i]);
                        new_dest.push_back(this->snode[this_neigh].sedges[this_mrg].dest_node[i]);
                        new_mults.push_back(this->snode[this_neigh].sedges[this_mrg].mults[i]);
                    }                
                }
                this->snode[this_neigh].sedges.erase(this_mrg);
            } 
        }
        if(is_pointed)
        {
            //cout << "Addig " << new_node << " to " << this_neigh <<"\n";
            this->snode[this_neigh].pointed_by.insert(new_node);
        } 
    }

#if PROFILE_MERGE
    gettimeofday( &stop, NULL);
    dur =  stop.tv_usec - start.tv_usec;
    cout << "\nphase 3: " << (size_t) dur;
#endif

    /*
    NOTE: At this point all super-edge connections have been redirected to the new node (old nodes have been cut off)
    */

    /**********************************************************
    PHASE 4
    ***********************************************************/
#if PROFILE_MERGE   
    gettimeofday( &start, NULL);
#endif    
    /*
        Remove merged supernodes
    */
    for(auto& merged : mrg_ids)
    {
        this->ids.erase(merged);
        this->snode.erase(merged);
    }
    this->nnodes -= mrg_ids.size();

#if PROFILE_MERGE
    gettimeofday( &stop, NULL);
    dur =  stop.tv_usec - start.tv_usec;
    cout << "\nphase 4: " << (size_t) dur;
#endif

    /**********************************************************
    PHASE 5
    ***********************************************************/
#if PROFILE_MERGE   
    gettimeofday( &start, NULL);
#endif    
    /*
    Update summary layer information: 1) Glyph, 
                                      2) Representative multiplicities,
                                      3) Active super-edges
                                      4) Costs  
    */

   /*
     Get representative multiplicities and their costs
   */
    this->snode[new_node].decide_glyph();
    this->snode[new_node].get_rewiring_cost();
    this->snode[new_node].get_rep_mult();

    for(auto it = this->snode[new_node].sedges.begin(); it != this->snode[new_node].sedges.end(); ++it)
    {
        it->second.get_rep_mult();
        it->second.get_rewiring_cost(this->snode[new_node].glyph, this->snode[it->first].glyph,
                                     this->snode[new_node].size, this->snode[it->first].size,
                                     this->snode[new_node].hub, this->snode[it->first].hub,
                                     this->snode[new_node].nodeset, this->snode[it->first].nodeset);        
    }    
    this->snode[new_node].update_edge_deltas();

#if PROFILE_MERGE
    gettimeofday( &stop, NULL);
    dur =  stop.tv_usec - start.tv_usec;
    cout << "\nphase 5a: " << (size_t) dur;  

    gettimeofday( &start, NULL);
#endif

    for(auto& this_neigh : neigh)
    {
        /*
            Multiplicities of merged_neighborhood superedges
        */
        if(this->snode[new_node].pointed_by.count(this_neigh) > 0)
        {
            for(auto it = this->snode[this_neigh].sedges.begin(); it != this->snode[this_neigh].sedges.end(); ++it)
            {
                if(it->first == new_node)
                {

                    it->second.get_rep_mult();
                    it->second.get_rewiring_cost(this->snode[this_neigh].glyph, this->snode[new_node].glyph,
                                                 this->snode[this_neigh].size, this->snode[new_node].size,
                                                 this->snode[this_neigh].hub, this->snode[new_node].hub,
                                                 this->snode[this_neigh].nodeset, this->snode[new_node].nodeset);   
                    // if(new_node == 256 && this_neigh == 22)
                    // {
                    //  cout << " HERE !\n";                              
                    // cout << " edge_delta " << it->second.rewiring_cost - it->second.non_rewiring_cost << "\n";                                                                          
                    // }
                }

            }
            this->snode[this_neigh].update_edge_deltas();
        }
    }     

#if PROFILE_MERGE
    gettimeofday( &stop, NULL);
    dur =  stop.tv_usec - start.tv_usec;
    cout << "\nphase 5b: " << (size_t) dur;  

    gettimeofday( &start, NULL);
#endif

    /*
        Update and sort edge_deltas 
    */
    this->snode[new_node].update_edge_deltas();
    this->snode[new_node].activate_edges(this->nnodes);
    for(auto it = this->snode[new_node].pointed_by.cbegin(); it != this->snode[new_node].pointed_by.cend(); ++it )
    {
        this->snode[*it].update_edge_deltas();
        this->snode[*it].activate_edges(this->nnodes);
    }    

#if PROFILE_MERGE
    gettimeofday( &stop, NULL);
    dur =  stop.tv_usec - start.tv_usec;
    cout << "\nphase 5c: " << (size_t) dur;  
#endif

//    gettimeofday( &start, NULL);

    /*
        Re-evaluate edge activation
        NOTE: This is to deal both with local changes in connectivity, as well as global shrinking of #nnodes
        UPDATE 11/26: Using the approximate edge activation rule, I only need to look at local nodes 
    */
    // for(auto it = this->snode.begin(); it != this->snode.end(); ++it)
    //     it->second.activate_edges(this->nnodes);

    // gettimeofday( &stop, NULL);
    // dur =  stop.tv_usec - start.tv_usec;
    // cout << "\nphase 5d: " << (size_t) dur;    

    return new_node;
}



/******************************************************************************************************************
 Expand supernodes to the simple nodes and edges they contain
*******************************************************************************************************************/
void graph::expand_snodes(vector<int>& exp_ids)
{
    if(exp_ids.empty())
        return;

    unordered_set<int> exp_set(exp_ids.cbegin(), exp_ids.cend());

    /* Find the set of snodes that points to the expanded snodes (excluding themselves) */
    unordered_set<int> neigh;
    for(auto& n : exp_ids)
        neigh.insert(this->snode[n].pointed_by.cbegin(), this->snode[n].pointed_by.cend());
    for(auto& n : exp_ids)
        neigh.erase(n);    


    /**********************************************************
    PHASE 1
    ***********************************************************/    

    /*
        Expand nodes to simple contained nodes
    */
    for(auto& n : exp_ids)
    {
        /* Enter a new snode for every node contained */
        this->ids.insert( this->snode[n].nodeset.cbegin(), this->snode[n].nodeset.cend());
        
        for(auto it = this->snode[n].nodeset.cbegin(); it != this->snode[n].nodeset.cend(); ++it)
        {
            this->snode.insert( pair<int,node>(*it, node()));
            this->snode[*it].type = this->snode[n].type;
            this->snode[*it].nodeset.insert(*it);
            this->nnodes++;
        }    

        /* Break down node sets and ingraphs */
        for(auto it = this->snode[n].ingraph.cbegin(); it != this->snode[n].ingraph.cend(); ++it)
        {
            int node_from = it->first; 
            for(auto it2 = it->second.cbegin(); it2 != it->second.cend(); ++it2)
            {
                int node_to = it2->first;
                size_t mult = it2->second;

                //cout << " ingraph " << node_from << " " << node_to <<"\n";

                if(node_from != node_to)
                {
                    /* Take ingraph edges and put them into the expanded super-edges */
                    this->snode[node_from].sedges.insert( pair<int,edge>( node_to, edge()));
                    this->snode[node_from].sedges[node_to].src_node.push_back(node_from);
                    this->snode[node_from].sedges[node_to].dest_node.push_back(node_to);
                    this->snode[node_from].sedges[node_to].mults.push_back(mult);

                    /* Add super-edge delta */
                    this->snode[node_from].sedge_deltas.reserve(this->snode[node_from].sedge_deltas.size() + 1);
                    this->snode[node_from].sedge_deltas.push_back(pair<int, double>(node_to, -1.0) );    
                        
                    /* Increase multiplicities and related initial costs */
                    this->snode[node_from].sedges[node_to].rep_mult = mult;
                    this->snode[node_from].sedges[node_to].rep_mult_cost = 2.0*log2( (double) this->snode[node_from].sedges[node_to].rep_mult)  +1.0;
                    this->snode[node_from].sedges[node_to].non_rep_mult_cost = this->snode[node_from].sedges[node_to].rep_mult_cost;

                    /* Also add the pointed-by sets.. */
                    this->snode[node_to].pointed_by.insert(node_from);
                }
                else
                {
                    /* This is a self loop so it goes to the (single-node) ingraph */     
                    this->snode[node_from].self_loop = true;
                    this->snode[node_from].ingraph.insert(pair<int,unordered_map<int,size_t>>(node_from, unordered_map<int,size_t>()));
                    this->snode[node_from].ingraph[node_from].insert(pair<int,size_t>(node_from, mult));
                    this->snode[node_from].sl_cost = 0.0;
                    this->snode[node_from].non_sl_cost = 1.0;
                    this->snode[node_from].sedge_deltas.reserve(this->snode[node_from].sedge_deltas.size() + 1);
                    this->snode[node_from].sedge_deltas.push_back(pair<int, double>(SELF_ID, -1.0) );
                    this->snode[node_from].sl_mult = mult;
                    this->snode[node_from].sl_mult_cost = 2.0*log2( (double) this->snode[node_from].sl_mult) + 1.0;
                    this->snode[node_from].non_sl_mult_cost = this->snode[node_from].sl_mult_cost;
             
                }
            
            }
        }    
        this->nnodes--;
    }

    /* Expand the super-edges */
    for(auto& n : exp_ids)
    {   
        for(auto it = this->snode[n].sedges.cbegin(); it != this->snode[n].sedges.cend(); ++it)
        {                 
            int sedge_points_to = it->first;
            auto& this_sedge = it->second;
            if(exp_set.count(sedge_points_to) > 0 )
            {
                /* These edges are pointing within the expanded set */
                for(size_t i=0; i< this_sedge.size; i++ )
                {
                    int node_from = this_sedge.src_node[i];
                    int node_to = this_sedge.dest_node[i];
                    size_t mult = this_sedge.mults[i];

                    //cout << " sedges_within " << node_from << " " << node_to <<"\n";

                    /* Add as simple superedge */    
                    this->snode[node_from].sedges.insert( pair<int, edge >( node_to, edge()));
                    this->snode[node_from].sedges[node_to].src_node.push_back(node_from);
                    this->snode[node_from].sedges[node_to].dest_node.push_back(node_to);
                    this->snode[node_from].sedges[node_to].mults.push_back(mult);

                    /* Add super-edge delta */
                    this->snode[node_from].sedge_deltas.reserve(this->snode[node_from].sedge_deltas.size() + 1);
                    this->snode[node_from].sedge_deltas.push_back(pair<int, double>(node_to, -1.0) );    
                        
                    /* Increase multiplicities and related initial costs */
                    this->snode[node_from].sedges[node_to].rep_mult = mult;
                    this->snode[node_from].sedges[node_to].rep_mult_cost = 2.0*log2( (double) this->snode[node_from].sedges[node_to].rep_mult)  +1.0;
                    this->snode[node_from].sedges[node_to].non_rep_mult_cost = this->snode[node_from].sedges[node_to].rep_mult_cost;                    

                    /* Also change the pointed-by sets of dest */
                    this->snode[node_to].pointed_by.insert(node_from);

                }
            }
            else
            {
                /* Remove the pointed-by reference from remote supernode */
                this->snode[sedge_points_to].pointed_by.erase(n);

                /* These edges point out of the expanded set  */
                for(size_t i=0; i< this_sedge.size; i++ )
                {
                    int node_from = this_sedge.src_node[i];
                    int node_to = this_sedge.dest_node[i];
                    size_t mult = this_sedge.mults[i];

                    //cout << " sedges_out " << node_from << " " << node_to <<"\n";

                    if(this->snode[node_from].sedges.count(sedge_points_to) == 0)
                    {
                        this->snode[node_from].sedges.insert(pair<int, edge>( sedge_points_to , edge() ));
                        this->snode[node_from].sedges[sedge_points_to].src_node.push_back(node_from);
                        this->snode[node_from].sedges[sedge_points_to].dest_node.push_back(node_to);
                        this->snode[node_from].sedges[sedge_points_to].mults.push_back(mult);

                        this->snode[node_from].sedge_deltas.reserve(this->snode[node_from].sedge_deltas.size() + 1);
                        this->snode[node_from].sedge_deltas.push_back(pair<int, double>(sedge_points_to, -1.0) );    

                        /* Increase multiplicities and related initial costs */
                        this->snode[node_from].sedges[sedge_points_to].rep_mult = mult;
                        this->snode[node_from].sedges[sedge_points_to].rep_mult_cost = 2.0*log2( (double) this->snode[node_from].sedges[sedge_points_to].rep_mult)  +1.0;
                        this->snode[node_from].sedges[sedge_points_to].non_rep_mult_cost = this->snode[node_from].sedges[sedge_points_to].rep_mult_cost;      
                    }
                    else
                    {
                        this->snode[node_from].sedges[sedge_points_to].simple_edge = false;
                        this->snode[node_from].sedges[sedge_points_to].size ++;

                        this->snode[node_from].sedges[sedge_points_to].src_node.reserve(this->snode[node_from].sedges[sedge_points_to].src_node.size() + 1);
                        this->snode[node_from].sedges[sedge_points_to].dest_node.reserve(this->snode[node_from].sedges[sedge_points_to].dest_node.size() + 1);
                        this->snode[node_from].sedges[sedge_points_to].mults.reserve(this->snode[node_from].sedges[sedge_points_to].mults.size() + 1);

                        this->snode[node_from].sedges[sedge_points_to].src_node.push_back(node_from);
                        this->snode[node_from].sedges[sedge_points_to].dest_node.push_back(node_to);
                        this->snode[node_from].sedges[sedge_points_to].mults.push_back(mult);
                    }
                         
                    /* Add to pointed-by references of remote supernode  */
                    this->snode[sedge_points_to].pointed_by.insert(node_from);
                
                }
            }

        }
    }

    /* Update costs and deltas for non-trivial outgoing super-edges */
    for(auto& n : exp_ids)
    {   
        for(auto& v : this->snode[n].nodeset)
        {
            for(auto it = this->snode[v].sedges.begin(); it != this->snode[v].sedges.end(); ++it)
            {       
                        it->second.get_rep_mult();
                        it->second.get_rewiring_cost(this->snode[v].glyph, this->snode[it->first].glyph,
                                                     this->snode[v].size, this->snode[it->first].size,
                                                     this->snode[v].hub, this->snode[it->first].hub,
                                                     this->snode[v].nodeset, this->snode[it->first].nodeset);  
            }
            this->snode[v].update_edge_deltas();
            this->snode[v].activate_edges(this->nnodes);
        }
    }            

    /**********************************************************
    PHASE 2
    ***********************************************************/    

    /*
        Expand incoming super-edges from neighbors
    */
    for(auto& n : neigh)
    {
        for(auto it = this->snode[n].sedges.begin();it != this->snode[n].sedges.end(); ++it)
        {
            if(exp_set.count(it->first))
            {
                /* This is incoming sedge to an expanded snode and needs to be broken */
                for(size_t i=0; i< it->second.size; i++)
                {

                    int node_from = it->second.src_node[i];
                    int node_to = it->second.dest_node[i];
                    size_t mult = it->second.mults[i];

                    if(this->snode[n].sedges.count(node_to) == 0)
                    {
                        this->snode[n].sedges.insert(pair<int, edge>( node_to , edge() ));
                        this->snode[n].sedges[node_to].src_node.push_back(node_from);
                        this->snode[n].sedges[node_to].dest_node.push_back(node_to);
                        this->snode[n].sedges[node_to].mults.push_back(mult);

                        this->snode[n].sedge_deltas.reserve(this->snode[node_from].sedge_deltas.size() + 1);
                        this->snode[n].sedge_deltas.push_back(pair<int, double>(node_to, -1.0) );    

                        /* Increase multiplicities and related initial costs */
                        this->snode[n].sedges[node_to].rep_mult = mult;
                        this->snode[n].sedges[node_to].rep_mult_cost = 2.0*log2( (double) this->snode[node_from].sedges[node_to].rep_mult)  +1.0;
                        this->snode[n].sedges[node_to].non_rep_mult_cost = this->snode[node_from].sedges[node_to].rep_mult_cost;      
                    }
                    else
                    {
                        this->snode[n].sedges[node_to].simple_edge = false;
                        this->snode[n].sedges[node_to].size ++;

                        this->snode[n].sedges[node_to].src_node.reserve(this->snode[n].sedges[node_to].src_node.size() + 1);
                        this->snode[n].sedges[node_to].dest_node.reserve(this->snode[n].sedges[node_to].dest_node.size() + 1);
                        this->snode[n].sedges[node_to].mults.reserve(this->snode[n].sedges[node_to].mults.size() + 1);

                        this->snode[n].sedges[node_to].src_node.push_back(node_from);
                        this->snode[n].sedges[node_to].dest_node.push_back(node_to);
                        this->snode[n].sedges[node_to].mults.push_back(mult);
                    }
                         
                    /* Add to pointed-by references of remote supernode  */
                    this->snode[node_to].pointed_by.insert(n);                   
                }
            }
        }
    }

    /* Erase original super-edges that point to expanded set */
    for(auto& n:neigh)
        for(auto& e: exp_set)
            this->snode[n].sedges.erase(e);

    /* Update costs and deltas for non-trivial outgoing super-edges */
    for(auto& v : neigh)
    {   
            for(auto it = this->snode[v].sedges.begin(); it != this->snode[v].sedges.end(); ++it)
            {       
                        it->second.get_rep_mult();
                        it->second.get_rewiring_cost(this->snode[v].glyph, this->snode[it->first].glyph,
                                                     this->snode[v].size, this->snode[it->first].size,
                                                     this->snode[v].hub, this->snode[it->first].hub,
                                                     this->snode[v].nodeset, this->snode[it->first].nodeset);  
            }
            this->snode[v].update_edge_deltas();
            this->snode[v].activate_edges(this->nnodes);
    } 


    /**********************************************************
    PHASE 3
    ***********************************************************/    
    /*
        Delete expanded super-nodes
    */
    for(auto& n : exp_ids)
        this->snode.erase(n);        

}



/*
    Obtain total encoding cost of summary
*/
pair<double,double> graph::get_cost(void)
{

    double scost = 0.0;

    /*
        SUMMARY COST
    */

    scost += 2.0*log2( (double) this->nnodes ) + 1.0;

    scost += 2.0*log2( (double) this->ntypes ) + 1.0;

    scost += (double) this->nnodes * log2( (double) this->ntypes );


    size_t nedges = 0;

    for(auto it = this->snode.cbegin(); it != this->snode.cend(); ++it)
    {
        scost += 2.0*log2( (double) it->second.size ) + 1.0;

        if(it->second.size > 1 )    
        {
            if(!this->no_mults)
                scost += 2.0*log2( (double) it->second.nd_mult) + 1.0;
        
            scost += 2.0;
        }

        size_t outdeg = 0;
        if( it->second.self_loop )
        {
            outdeg ++;
            if(!this->no_mults)
                scost += 2.0*log2( (double) it->second.sl_mult) + 1.0;
        }    

        for(auto it2 = it->second.sedges.cbegin(); it2 != it->second.sedges.cend(); ++it2)    
            if(it2->second.active)
            {
                outdeg++;
                if(!this->no_mults)
                    scost += 2.0*log2( (double) it2->second.rep_mult ) + 1.0;
            }    

        scost += binom_cost( this->nnodes, outdeg ) + log2(this->nnodes) - 1.0 - ((outdeg>0)? 2.0*log2(this->nnodes) : 0.0); //<--- CHANGED
        nedges += outdeg;
    }

    /*
        CORRECTIONS COST
    */

    double rcost = 0.0;

    for(auto it = this->snode.cbegin(); it != this->snode.cend(); ++it)
    {

        /*
            Glyph cost
        */

        rcost += it->second.glyph_cost;
        if(!this->no_mults)
            rcost += it->second.nd_mult_cost;
        // cout << "Glyph cost: " << it->second.glyph_cost <<"\n\n";

        /*
            Align cost
        */
        rcost += binom_cost( this->nnodes_full, it->second.size ) - 2.0*log2(it->second.size) - 1.0;
        if(it->second.glyph == INSTAR || it->second.glyph == INSTAR)
            rcost += log2( (double) it->second.size);


        /*
            Self loop correction cost
        */
        if( it->second.self_loop )
        {
            rcost += it->second.sl_cost;
            if(!this->no_mults)
                rcost += it->second.sl_mult_cost;
        }           
        else
        {
            rcost += it->second.non_sl_cost;
            if(!this->no_mults)
                rcost += it->second.non_sl_mult_cost;            
        }
        

        /*
            For every adjacent edge/superedge active/non-active
        */
        for(auto it2 = it->second.sedges.cbegin(); it2 != it->second.sedges.cend(); ++it2)    
            if(it2->second.active)
            {
                rcost += it2->second.rewiring_cost;
                if(!this->no_mults)
                    rcost += it2->second.rep_mult_cost;
            }    
            else
            {
                rcost += it2->second.non_rewiring_cost;
                if(!this->no_mults)
                    rcost += it2->second.non_rep_mult_cost;
            }
            
    }


    return pair<double,double>(scost, rcost);
}




/*  Obtain  simple copy of forward super-graph topology  */
unordered_map<int,vector<int>> graph::get_fwd_graph(void)
{
    unordered_map<int,vector<int>> res;

    for(auto fromit = this->snode.cbegin(); fromit != this->snode.cend(); ++fromit)
    {
        int from = fromit->first;
        if(res.count(from) == 0)
            res.insert( pair<int,vector<int>>( from, vector<int>()));

        res[from].reserve(fromit->second.sedges.size()+1);

        for(auto toit = fromit->second.sedges.cbegin(); toit != fromit->second.sedges.cend(); ++toit)
            res[from].push_back(toit->first);

#if(EXPORT_SELF_LOOPS)
        if(fromit->second.self_loop)
            res[from].push_back(from);
#endif

        sort(res[from].begin(), res[from].end()); 
    }

    assert(res.size() == (size_t) this->nnodes );

    return res;
}


/*  Obtain  simple copy of inverse super-graph topology  */
unordered_map<int,vector<int>> graph::get_rev_graph(void)
{

    unordered_map<int,vector<int>> res;

    for(auto toit = this->snode.cbegin(); toit != this->snode.cend(); ++toit)
    {
        int to = toit->first;
        if(res.count(to) == 0)
            res.insert( pair<int,vector<int>>( to, vector<int>(toit->second.pointed_by.cbegin(),toit->second.pointed_by.cend())));

#if(EXPORT_SELF_LOOPS)
        if(toit->second.self_loop)
            res[to].push_back(to);
#endif
        sort(res[to].begin(), res[to].end());    
    }

    assert(res.size() == (size_t) this->nnodes );

    return res;

}



/* 
    Cleanup empty supernodes 
*/
void graph::cleanup(void)
{
    set<int> empty_snodes;
    for(auto it= this->snode.begin(); it != this->snode.end(); ++it)
        if(it->second.nodeset.empty())
            empty_snodes.insert(it->first);
            
    for(auto& v: empty_snodes)
    {
        this->snode.erase(v);
        this->ids.erase(v);
    }

    this->nnodes = this->snode.size();
}




/* 
    Save summary graph to file
*/
void graph::save_as(string fname)
{

    unordered_map<int, string> glyph_map;
    glyph_map.insert( pair<int, string>( (int) CLIQUE, "CLIQUE"));
    glyph_map.insert( pair<int, string>( (int) DISC, "DISC"));    
    glyph_map.insert( pair<int, string>( (int) INSTAR, "INSTAR"));    
    glyph_map.insert( pair<int, string>( (int) OUTSTAR, "OUTSTAR"));
    glyph_map.insert( pair<int, string>( (int) NONE, "NONE"));

    ofstream myfile;

    myfile.open(fname);

    myfile << "# NODES\n";
    myfile << "# ID  Size  Type  Glyph Glyph-rep-mult Self-loop SL-rep-mult\n";
    for(auto it = this->snode.cbegin(); it != this->snode.cend(); ++it)
    {
        myfile << it->first << " ";
        myfile << it->second.size << " ";
        myfile << it->second.type << " ";
        myfile << glyph_map[(int) it->second.glyph] << " ";
        myfile << it->second.nd_mult << " ";
        myfile << it->second.self_loop << " ";
        if(it->second.self_loop)
            myfile << it->second.sl_mult << "\n";
        else
        {
            myfile << "0\n";
        }
    }

    myfile << "\n";

    myfile << "# EDGES\n";
    myfile << "# Source-ID  Dest-ID Rep-mult\n";
    for(auto it = this->snode.cbegin(); it != this->snode.cend(); ++it)
    {
        for(auto it2 = it->second.sedges.cbegin(); it2 != it->second.sedges.cend(); it2++ )
        {
            if(it2->second.active)
            {
                myfile << it->first << " ";
                myfile << it2->first << " ";
                //cout << it2->second.size << " ";
                myfile << it2->second.rep_mult << "\n";
            }
        }
    }

    myfile.close();
}


