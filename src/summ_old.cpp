#include <cstdlib>
#include <iostream>
#include <ostream>
#include <iterator>
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <random>
#include <iomanip>
#include <unordered_set>
#include <unordered_map>
#include <cmath>
#include <boost/functional/hash.hpp>

#include "omp.h"

#include "utils.h"
#include "summ.h"
#include "graph.h"

#define CHECK_FP false

#define GET_SIM true

#define EPS 100.0

#define STAR_THRESHOLD_1 0.8

#define STAR_THRESHOLD_2 0.2

using namespace std;

typedef mt19937 rng_type;



/*------------------------------------------------------------------------

                        SUMMARIZATION ALGORITHM(S)

-------------------------------------------------------------------------*/

/*
    Generate r min-hash signatures for every node.
    Then hash the signatures.
    Return hash table
*/
static void generate_hash_table(unordered_map<string, set<int>>& htable , unordered_map<int,vector<int>>& adj, size_t r)
{

    /* Random number generator from uniform distribution */
    uniform_int_distribution<rng_type::result_type> udist(0, adj.size() -1);
    rng_type rng;    

    vector<int> temp(r, 0);

    /*
        TODO: May do this more efficiently with OMP but will need an array since res.insert() probably NOT thread safe
    */

    /* Generate a random offset of deterministic seed values
       This is to make this hash-table different than the previous one */
    int rand_offset = rand();

    /* Iterate over nodes */
    for(auto it = adj.cbegin(); it != adj.cend(); ++it)
    {
        /* Initialize min-hash signatures to maximum value */
        fill(temp.begin(), temp.end(), adj.size() );

        /* Generate signatures for this node */
        for(auto& v : it->second )
        {
            /* Seed random with current item value */
            rng.seed( (rng_type::result_type const) (v + rand_offset) );
            
            for(size_t i =0; i< temp.size(); i++)
                temp[i] = min(temp[i], (int) udist(rng));
        }


        stringstream key;
        copy(temp.cbegin(), temp.cend(), ostream_iterator<int>(key, ""));    
        
         if( htable.count(key.str()) == 0)
             htable.insert( pair<string, set<int>>( key.str(), set<int>()));
        htable[key.str()].insert( it->first); 
    }
      
}



/*
    Hash (LSH) nodes into sets (sorted vectors).
    Return ranked list of sets according to their scores.
    Each entry is a pair of set and score of set  
    Parameter: size of band r, number of bands b
*/
vector<pair< pair<size_t,double> ,set<int>>> increment_lsh(unordered_map<int,vector<int>>& adj, size_t r, size_t b )
{

    vector<pair<pair<size_t,double>,set<int>>> res;
    res.reserve(adj.size()*b);  

    /* This maps node ids to cluster ids */
    unordered_map<int, int> nid2cid;

    /* This is maps cluster id to node ids contained */
    unordered_map<int, set<int>> cid2nid;

    /* Cluster ids initialized 1:N*/
    int cid = 1;
    for(auto it = adj.cbegin(); it != adj.cend(); ++it)
        nid2cid.insert(pair<int, int>(it->first, cid++));

    /* Seed here once*/
    srand( 1 );

    unordered_map<string, set<int>> htable;

    for(size_t i = 0; i< b; i++)
    {
        /* This is the approximate value of the threshold accod=rding to LSH theory */
        double threshold =  pow(1.0/(double) ( 1 + i ), 1.0/(double)r );

        /* Populate hash table using r min-hashes */
        generate_hash_table( htable, adj, r);

#if(false)
        cout << "\n\nHtable\n\n";
        for(auto it = htable.cbegin(); it != htable.cend(); ++it )
        {
            cout << "\nBucket " << it->first << ": " ;
            for(auto it2 = it->second.cbegin(); it2 != it->second.cend(); ++ it2)
                cout << " " << *it2;    
            cout << "\n";    
        }       
#endif

        /* Iterate over buckets with more than one node */
        for(auto it = htable.cbegin(); it != htable.cend(); ++it )
        {
            if(it->second.size()>1)
            {   
#if(CHECK_FP)                
                /* Check for false positive pairs */
                cout << "\n\n # FP-ratio: " << check_FP(it->second, adj, threshold) << "\n";
#endif

                /* Update node to cluster mapping */
                int min_cid = adj.size();
                for(auto it2 = it->second.cbegin(); it2 != it->second.cend(); ++ it2)
                    min_cid = ( nid2cid[*it2] < min_cid ) ? nid2cid[*it2] : min_cid;
                for(auto it2 = it->second.cbegin(); it2 != it->second.cend(); ++ it2)
                    nid2cid[*it2] = min_cid;

                /* Add new cluster OR merge one of the existing ones 
                    This to update our clustering "working-space" */    
                if(cid2nid.count(min_cid) == 0)
                    cid2nid.insert( pair<int, set<int>>( min_cid, set<int>() ));

                cid2nid[min_cid].insert(it->second.cbegin(), it->second.cend());                

                /* Check if an identical to this set has be stored before*/
                bool identical_set = false;
                for(auto setit = res.cbegin(); setit != res.cend(); ++setit )
                {
                    if(cid2nid[min_cid].size() == setit->second.size())
                        if(cid2nid[min_cid] == setit->second)
                        {
                            identical_set = true;
                            break;
                        }
                            
                }

                if(!identical_set)
                {
#if(GET_SIM)
                    /* Get true average similarity of nodes in cluster */
                    double sim_ = get_sim(cid2nid[min_cid], adj); 
#else
                    /* Just use the current threshold ignoring the false positives */
                    double sim_ = threshold;                    
#endif
                    /* Store new cluster (id is not needed) to results 
                    Score will be added yet */
                    res.push_back( pair<pair<size_t,double>,set<int>>( pair<size_t, double>(cid2nid[min_cid].size(), sim_ ), set<int>(cid2nid[min_cid].cbegin(),cid2nid[min_cid].cend())));

#if(true)
                    cout << "\n\nAdded cluster (size: " << res.back().first.first 
                        << ", sim: " << res.back().first.second << ", exp " << (double)(r+1)/(double)(r+2) <<" )\n\n";
                    for(auto& v : res.back().second )
                        cout << " " << v;    
                    cout << "\n";    
#endif                
                }
            }    
        }      

#if(false)
        cout << "\n\nCluster mapping\n\n";
        cout << "\nnid: " ;
        for(auto it = nid2cid.cbegin(); it != nid2cid.cend(); ++ it)
            cout << " " << it->first;    
        cout << "\ncid: " ;
        for(auto it = nid2cid.cbegin(); it != nid2cid.cend(); ++ it)
            cout << " " << it->second;    
        cout << "\n";    

#endif

        /* Clear hash table for new iteration */
        htable.clear();
    }


    res.shrink_to_fit();
    return res;
}



/*
Merges forward and reverse lists.
Sorts sets according to size/similarity.
Returns top-K candidate sets
*/
vector<set<int>> get_candidates(vector<pair<pair<size_t, double>,set<int>>>& fwd_list,
                                          vector<pair<pair<size_t, double>,set<int>>>& rev_list,
                                          unordered_map<int,vector<int>>& adj_fwd,
                                          unordered_map<int,vector<int>>& adj_rev,                                          
                                          size_t K)
{

    set<int> sinks;
    set<int> sources;
    set<int> disc;

    /* Put aside sink nodes (do not point anywhere) */
    for(size_t i=0; i< fwd_list.size(); i++ )
        if(isnan(fwd_list[i].first.second))
        {
            sinks.insert(fwd_list[i].second.cbegin(), fwd_list[i].second.cend());
            fwd_list.erase( fwd_list.begin() + i);
        }
            
     /* Put aside source nodes (are not pointed) */       
    for(size_t i=0; i< rev_list.size(); i++ )
        if(isnan(rev_list[i].first.second))
        {
            sources.insert(rev_list[i].second.cbegin(), rev_list[i].second.cend());
            rev_list.erase( rev_list.begin() + i);
        }

    /* Find completely incative nodes (neither point, nor pointed, SELF-LOOPs possible)  */        
    set_intersection(sources.begin(), sources.end(), sinks.begin(), sinks.end(), inserter( disc, disc.begin()));

    /* Clean sinks and sources from completely disconected nodes */
    for(auto& v : disc)
    {
        sinks.erase(v);
        sources.erase(v);
    }

    /* Check for potential INSTAR spokes
       This will be nodes with similar outgoing connections and little-to-no incoming connections
       If they have 1 common destination, add it to the set (since it will be the hub)*/
    set<int> temp;   
    for(auto it = fwd_list.begin(); it != fwd_list.end(); ++it)
    {
        /* Check if outgoing similarity is above a threhold */
        if(it->first.second > STAR_THRESHOLD_1)
        {
            set_intersection(it->second.begin(), it->second.end(), sources.begin(),
                             sources.end(), inserter(temp, temp.begin()));
            /* Check if proposrtion of sources is large enough */
            if( (double) temp.size() / (double) it->second.size() > STAR_THRESHOLD_2 )
            {
                /* OK this is probably an INSTAR..
                   Find the hub */
                int hub;
                auto last_intersection = adj_fwd[*it->second.begin()];
                vector<int> current_intersection;
                auto it2 = it->second.begin();
                it2++;
                for(; it2 != it->second.end(); ++it2 )
                {
                    set_intersection(last_intersection.begin(), last_intersection.end(),
                                     adj_fwd[*it2].begin(), adj_fwd[*it2].end(),
                                     back_inserter(current_intersection));
                    swap(last_intersection, current_intersection);
                    current_intersection.clear();
                }

                /* If hub is found, insert it */
                if(last_intersection.size() == 1)
                {
                    hub = last_intersection[0];
                    it->second.insert(hub);
                    it->first.first++;
                    // cout << "\n\n Found INSTAR with hub: " << hub
                    //      << ", size: " << it->second.size() << ", score " << it->first.second << "\n\n"; 
                } 
            }
        }    
    }

    /* Similar for OUTSTAR */
    for(auto it = rev_list.begin(); it != rev_list.end(); ++it)
    {
        /* Check if incoming similarity is above a threhold */
        if(it->first.second > STAR_THRESHOLD_1)
        {
            set_intersection(it->second.begin(), it->second.end(), sources.begin(),
                             sinks.end(), inserter(temp, temp.begin()));
            /* Check if proposrtion of sinks is large enough */
            if( (double) temp.size() / (double) it->second.size() > STAR_THRESHOLD_2 )
            {
                /* OK this is probably an OUTSTAR..
                   Find the hub */
                int hub;
                auto last_intersection = adj_rev[*it->second.begin()];
                vector<int> current_intersection;
                auto it2 = it->second.begin();
                it2++;
                for(; it2 != it->second.end(); ++it2 )
                {
                    set_intersection(last_intersection.begin(), last_intersection.end(),
                                     adj_rev[*it2].begin(), adj_rev[*it2].end(),
                                     back_inserter(current_intersection));
                    swap(last_intersection, current_intersection);
                    current_intersection.clear();

                    /* If hub was found, insert it */
                    if(last_intersection.size() == 1)
                    {
                        hub = last_intersection[0];
                        it->second.insert(hub);
                        it->first.first++;
                        // cout << "\n\n Found OUTSTAR with hub: " << hub
                        //      << ", size: " << it->second.size() << ", score " << it->first.second << "\n\n";                         
                    } 
                }
            }
        }                 
    }

    /* Merge lists */
    vector<set<int>> res;
    res.push_back(disc);
    res.reserve( fwd_list.size() + rev_list.size());
    for(auto it = fwd_list.cbegin();it != fwd_list.cend(); ++it )
        res.push_back( it->second );
    for(auto it = rev_list.cbegin();it != rev_list.cend(); ++it )
        res.push_back( it->second );        

    /* Sort */

    /* Add disconnected set in the beginning since it will almost definetly help reduce cost */


    /* Return candidates */
    return res;
}



/*
 Measure cost reduction from candidate set merging 
 */
double cost_reduction( double cost_before, graph& graph_,  set<int>& cand_set )
{
    if(cand_set.size() <= 1)
        return 0;

    /* Separate candidate set if types dont agree */
    unordered_map< int, vector<int> > cand_per_type;
    for(auto& v : cand_set)
    {
        if(cand_per_type.count(graph_.snode[v].type) == 0 )
            cand_per_type.insert( pair<int, vector<int>>( graph_.snode[v].type, vector<int>()));
        cand_per_type[graph_.snode[v].type].reserve(cand_per_type[graph_.snode[v].type].size() + 1 );
        cand_per_type[graph_.snode[v].type].push_back(v);
    }
        
    /* Merge sets and save ids of new super-nodes */
    vector<int> new_supernodes;
    new_supernodes.reserve(cand_per_type.size());
    for(auto it = cand_per_type.begin(); it != cand_per_type.end(); ++it)
        if(it->second.size() > 1)
            new_supernodes.push_back(graph_.merge_snodes( it->second ));

    graph_.cleanup();

    /* Compute cost after merging */
    pair<double,double> costs_after = graph_.get_cost();
    double cost_after = costs_after.first + costs_after.second;

    /* Undo merging by expanding new super-nodes */
    graph_.expand_snodes(new_supernodes);    

    pair<double,double> c = graph_.get_cost();
    double ct = c.first + c.second;
    
    cout << "\nbefore " << cost_before << " rec " << ct << "\n";
    assert( fabs(ct - cost_before) <= EPS * fabs(ct) );

    return cost_before - cost_after;
}



/*
    Greedy summarization of graph given list of candidate sets.
    Summarizes graph in place.
    Returns total cost before and after summarization.
*/
pair<double, double> greedy_summ(graph& graph_,  vector<set<int>>& candidates)
{
    pair<double,double> res;

    unordered_map< int, vector<int> > cand_per_type;

    size_t iter_ = 0;

    while(true)
    {
        cout << "\n\nIteration: " << ++iter_;

        /* Get cost of graph at this iteration */
        pair<double,double> costs_before = graph_.get_cost();
        double cost_before = costs_before.first + costs_before.second;

        if(iter_ == 1)
            res.first = cost_before;

        /* Iterate over list of candidate sets to find best */
        size_t best_cand, cand_counter = 0;
        double best_delta = -1.0;
        for(auto& cand_set : candidates)
        {
            // cout << "\nCand set: ";
            // for(auto& v: cand_set)
            //     cout <<" " << v;
            // cout << "\n";

            double delta = cost_reduction( cost_before, graph_, cand_set );
            cout << " --> delta: " <<  delta <<"\n";
            if(delta > best_delta)
            {
                best_delta = delta;
                best_cand = cand_counter;
            }
            cand_counter++;
        }

        /* Terminate if no improvement possible */
        if( best_delta <= 0.0 )
        {
            cout << "\n\nNo further improvement possible..\n\n";
            break;
        }
            

        /* MERGE BEST CANDIDATE SET */
        ///////////////////////////////////////////////////////////
        cout << "\nBest Candidate: ";
        for(auto& v : candidates[best_cand]) cout << v << " ";
        cout << "\n";
        /* Separate candidate set if types dont agree */
        for(auto& v : candidates[best_cand])
        {
            if(cand_per_type.count(graph_.snode[v].type) == 0 )
                cand_per_type.insert( pair<int, vector<int>>( graph_.snode[v].type, vector<int>()));
            cand_per_type[graph_.snode[v].type].reserve(cand_per_type[graph_.snode[v].type].size() + 1 );
            cand_per_type[graph_.snode[v].type].push_back(v);
        }
            
        /* Merge sets per type */
        for(auto it = cand_per_type.begin(); it != cand_per_type.end(); ++it)
            if(it->second.size() > 1)
                graph_.merge_snodes( it->second );
        cand_per_type.clear();        
        /////////////////////////////////////////////////////////////

        //graph_.display();

        /* Subtract nodes that were summarized from other candidate sets */
        for(size_t i=0; i<candidates.size(); i++)
        {
            if(i!= best_cand)
            {
              for(auto& v : candidates[best_cand])   
                  candidates[i].erase(v); 
            }
        }


        /* Remove the candidate set that was used, and empty/singleton sets */
        candidates[best_cand].clear();
        for( auto it = candidates.begin(); it != candidates.end(); )
            if( it->size() <= 1)
                candidates.erase(it);
            else
                ++it;
                


        /* Also terminate if set becomes empty */
        if( candidates.size() == 0 )
        {
            cout << "\n\nExhausted candidate sets..\n\n";
            break;
        }
    }

    /* Get cost of graph after summarization */
    pair<double,double> costs_after = graph_.get_cost();
    double cost_after = costs_after.first + costs_after.second;

    res.second = cost_after;

    return res;

}