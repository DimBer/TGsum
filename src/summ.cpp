#include <cstdlib>
#include <iostream>
#include <ostream>
#include <iterator>
#include <vector>
#include <deque>
#include <memory>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <random>
#include <iomanip>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <cmath>
#include <boost/functional/hash.hpp>
#include <chrono>
#include <sys/time.h>

#include "omp.h"

#include "utils.h"
#include "summ.h"
#include "graph.h"

#define CHECK_FP false

#define GET_SIM true

#define EPS 100.0

#define STAR_THRESHOLD_1 0.25

#define STAR_THRESHOLD_2 0.01

#define CHECK_BUCKETS true

#define BUCKETIZE_SELF_LOOPS true

#define BUCKETIZE_EXCHANGE_PAIRS true

#define SAVE_CAND_GRAPHS true

#define PROFILE_MRG_EXP_SUB false

#define ONE_CLIQUE false

using namespace std;

typedef mt19937 rng_type;


/* 
  This will store all necessairy information on a node cluster
  1) dnode: Nodes contained sorted by outdegree + indegree (i.e., size of corresponding fwdrev list)
  2) hpair: Latent ``hiden'' pairs sorted by similarity
  3) Lower-bound on size of maximum clique contained
  Each cluster will be added to map and indexed by label value
*/
typedef struct cluster{
    deque<pair<size_t,int>> dnode;
    deque<pair<double,pair<int,int>>> hpair;
#if ONE_CLIQUE
    size_t lmaxq;
#endif
} Cluster;

/* 
 Jaccard coefficient between out-edges of two nodes
*/
inline double jaccard(unordered_map<int,vector<int>>& adj, int node_a, int node_b)
{
    double intersection_ = (double) intersection_size(adj[node_a], adj[node_b]);

    double union_ = (double) union_size(adj[node_a], adj[node_b]);

    if(intersection_ == 0.0 && union_ == 0.0 )
        return 0.0;

    return intersection_/union_;
}

/*
    Recursive HEURISTIC for approximating maximum clique on a subset of nodes U of undirected graph G.
    Clique is returned in Q. We may only look for cliques > lmaxq, so terminate early if such clique is impossible.
*/
inline void clique_heuristic(unordered_map<int, set<int>>& G, set<int>& U,
                             deque<int>& Q, size_t lmaxq)
{
    /* Pick new candiate clique node while U is not empty */
    while(!U.empty())
    {
        /* Terminate early if lower bound cant be reached */
        if( Q.size() + U.size() <= lmaxq )
            return;

        /* Pick u: maximum degree candidate from U */
        int u; size_t maxd=0;
        for(auto z : U)
        {
            size_t d_z = G[z].size();
            if(d_z > maxd)
            {
                maxd = d_z;
                u = z;
            } 
        } 

        /* Put u it into Q */
        Q.push_back(u);

        /* Remove from U nodes not adjacent to u 
           using intersection in place.. */    
        auto it1 = U.begin();
        auto it2 = G[u].begin();
        while ( (it1 != U.end()) && (it2 != G[u].end()) ) {
            if (*it1 < *it2) {
                U.erase(it1++);
            } else if (*it2 < *it1) {
                ++it2;
            } else { // *it1 == *it2
                    ++it1;
                    ++it2;
            }
        }

        /* Anything left in U from here on did not appear in G[u],
         so we remove it */
        U.erase(it1, U.end());
    }
}                             

/*
  Returns inverted index of candidate sets (i.e., maps every node to list of candidate sets it belongs)
*/
static unordered_map<int, deque<int>> cand_set_invert_index(deque<pair<pair<size_t, double>,set<int>>>& candidates)
{

    unordered_map<int, deque<int>> res;
    for(size_t i=0; i<candidates.size(); i++)
    {
        for(auto v: candidates[i].second )
        {
            if( res.count(v) == 0 )
                res.insert(pair<int,deque<int>>(v, deque<int>()));
            res[v].push_back(i);    
        }
    }
    return res;
}


/*
    Returns overlap candidate-set graph (i.e., two nodes have edge iff they correspoding sets overlap)
*/
static unordered_map<int, unordered_set<int>> get_cand_overlap_graph(unordered_map<int, deque<int>>& inverted_ind, size_t Nsets)
{
    unordered_map<int, unordered_set<int>> res;
    for(size_t i=0; i< Nsets; i++)
        res.insert( pair<int,unordered_set<int>>( i, unordered_set<int>()));

    for(auto it = inverted_ind.cbegin(); it!= inverted_ind.cend(); ++it)
    {
        size_t nsets = it->second.size();
        if(nsets > 1)
        {
            for(size_t i=0; i<nsets; i++)
            {
                for(size_t j=i+1; j<nsets; j++)
                {
                    res[it->second[i]].insert(it->second[j]);
                    res[it->second[j]].insert(it->second[i]);                    
                }
            }
        }
    }

    return res;
}


/*
    Returns proximity candidate-set graph (i.e., two nodes have edge iff they correspoding sets share an edge)
*/
static unordered_map<int, unordered_set<int>> get_cand_prox_graph(unordered_map<int, deque<int>>& inverted_ind,
                                                         unordered_map<int,vector<int>>& adj, size_t Nsets)
{

    unordered_map<int, unordered_set<int>> res;
    for(size_t i=0; i< Nsets; i++)
        res.insert( pair<int,unordered_set<int>>( i, unordered_set<int>()));

    for(auto it = adj.cbegin(); it!= adj.cend(); ++it)
    {
        int node_a = it->first;
        for(auto node_b : it->second)
        {
            if(node_a != node_b)
            {
                for(auto set_a : inverted_ind[node_a])
                {
                    for(auto set_b : inverted_ind[node_b])
                    {
                        if(set_a != set_b)
                        {
                            res[set_a].insert(set_b);
                            res[set_b].insert(set_a);                                    
                        }
                    }
                }
            }
        }
    }

    return res;
} 

/*
    Remove node from candidate graph
*/
static void cand_graph_remove_node(unordered_map<int, unordered_set<int>>& cand_graph, int node_ )
{
    if(cand_graph.count(node_) > 0)
    {
        for(auto n : cand_graph[node_])
        {
            cand_graph[n].erase(node_);
        }
        cand_graph.erase(node_);
    }
}


/*
    Kako-Haldroson (https://www.ru.is/~mmh/papers/WIS_WG.pdf) greedy
    approximation of the maximum weighted independent set
 */
static unordered_set<size_t> approx_wis(vector<double>& w, unordered_map<int, unordered_set<int>>& graph_, double& W)
{
    unordered_set<size_t> res;
    W = 0.0;

    /* Compute weighted degrees */
    vector<pair<size_t,double>> d_w(w.size());
    for(size_t i=0; i<w.size(); i++)
    {
        d_w[i].first = i;
        d_w[i].second = 0.0;
        for(auto v: graph_[i])
            d_w[i].second += w[v];
        d_w[i].second /= w[i];    
    }

    /* Sort nodes wrt weighted degrees */
    sort(d_w.begin(), d_w.end(), [](auto &left, auto &right) {return left.second < right.second;});

    /* Run greedy removal algorithm */
    size_t i=0;
    while(true)
    {
        cout << "\n\n check " << i << " size " << d_w.size();
        /* Find next node in graph with smallest weighted degree */
        while( graph_.count(d_w[i].first) == 0 && i < graph_.size())
            i++;

        /* Terminate if graph is empty */    
        if(i >= graph_.size())
            break;
        
        /* Put node in the wis */
        res.insert(d_w[i].first);    
        W += w[d_w[i].first];

        /* Remove node and its neighborhood */
        for(auto neigh: graph_[d_w[i].first])
        {
            for(auto neigh2: graph_[neigh])
                if(graph_.count(neigh2) > 1)
                    graph_[neigh2].erase(neigh);
            graph_.erase(neigh);    
        }
        graph_.erase(d_w[i].first);

        if(i++ >= graph_.size())
            break;
    }

    return res;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////

/*------------------------------------------------------------------------

                        CANDIDATE SET GENERATION

-------------------------------------------------------------------------*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
    Generate r min-hash signatures for every node.
    Then hash the signatures.
    Return hash table
*/
static void generate_hash_table(unordered_map<string, deque<int>>& htable , unordered_map<int,vector<int>>& adj, size_t r)
{

    /* Random number generator from uniform distribution */
    int range = INT32_MAX;
    uniform_int_distribution<rng_type::result_type> udist(0, range - 1 );
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
        /* Avoid hashing sinks/sources */
        if(!it->second.empty())
        {
            /* Initialize min-hash signatures to maximum value */
            fill(temp.begin(), temp.end(), range );

            /* Generate signatures for this node */
            for(auto& v : it->second )
            {
                /* Seed random with current item value */
                rng.seed( (rng_type::result_type const) (v + rand_offset) );
                
                for(size_t i =0; i< temp.size(); i++)
                    temp[i] = min(temp[i], (int) udist(rng));
            }

            /* Add node to hash table */
            stringstream key;
            copy(temp.cbegin(), temp.cend(), ostream_iterator<int>(key, ""));    

            // cout << "\nNode: " << it->first << " Key: " << key.str();
            
            if( htable.count(key.str()) == 0)
                htable.insert( pair<string, deque<int>>( key.str(), deque<int>()));
            htable[key.str()].push_back( it->first); 
        }
    }      
}



/*
    Hash (LSH) nodes into sets (sorted vectors).
    Return ranked list of sets according to their scores.
    Each entry is a pair of set and score of set  
    Parameter: size of band r, number of bands b
*/
deque<pair< pair<size_t,double> ,set<int>>> cand_gen(unordered_map<int,vector<int>>& adj, size_t r, size_t b_max )
{

    //////////////////////////////////////////////////////////////////////////////////////
    /* Initialize Clusters -- Label Array -- Similarity Graph */
    //////////////////////////////////////////////////////////////////////////////////////
    deque<pair<pair<size_t,double>,set<int>>> res;  

    /* Initialize node cluster labels (and clique sizes to 1) */
    unordered_map<int, pair<int,int>> nid2cid;
#if !ONE_CLIQUE    
    unordered_map<int, int> node_maxq;
#endif

    for(auto it = adj.cbegin(); it != adj.cend(); ++it)
    {
        nid2cid.insert(pair<int, pair<int,int>>(it->first, pair<int,int>(it->first, it->first)));
#if !ONE_CLIQUE    
        node_maxq.insert( pair<int, size_t>( it->first, 1 ));
#endif
    }

    /* Initialize clusters
       Map: Cluster IDs to clusters */
    unordered_map<int, cluster> cid2clust;
    for(auto it = adj.cbegin(); it != adj.cend(); ++it)
    {
        cluster temp;
        temp.dnode.push_back(pair<size_t,int>( it->second.size(), it->first));
        cid2clust.insert(pair<int, cluster>(it->first, temp));
#if ONE_CLIQUE
        temp.lmaxq = 1;
#endif
    }

    /* Initialize similarity graph */
    unordered_map<int, set<int>> gsim;
    for(auto it = adj.cbegin(); it != adj.cend(); ++it)
        gsim.insert( pair<int,set<int>>( it->first, set<int>()));

    /* 
       Initialize bucket map 
       (this is used to avoid checkign pairs that were in the same bucket)
    */
    unordered_map<int, int> nid2bucket;
    for(auto it = adj.cbegin(); it != adj.cend(); ++it)
        nid2bucket.insert(pair<int, int>(it->first, 0));

    //////////////////////////////////////////////////////////////////////////////////////

    /* MAIN LOOP */

    /////////////////////////////////////////////////////////////////////////////////////

    /* Seed here once*/
//    srand( 1 );

    unordered_map<string, deque<int>> htable;
    htable.rehash(1000000);

    /* Minimum threshold that we will consider overall */
    double min_thres = pow(1.0/(double) ( 1 + b_max ), 1.0/(double)r );

    for(size_t b = 0; b< b_max; b++)
    {
        /* This is the approximate value of the threshold accod=rding to LSH theory */
        double threshold =  pow(1.0/(double) ( 1 + b ), 1.0/(double)r );

        //cout << "\n\n\nThreshold: " << threshold << "\n";
        cout << "\n\n b = "<< b;

        /**********************************************************
         PHASE 1: Pop and add possible hidden pairs that satisfy new threshold
        ***********************************************************/   
        for(auto it= cid2clust.begin(); it!= cid2clust.end(); ++it)
        {
            size_t npairs = 0;
            while(!it->second.hpair.empty() && it->second.hpair.back().first > threshold)
            {
                int node_a = it->second.hpair.back().second.first;
                int node_b = it->second.hpair.back().second.second;
                gsim[node_a].insert(node_b);
                gsim[node_b].insert(node_a);
                it->second.hpair.pop_back();
                npairs++;
            }
            //if(npairs) cout << "\npoped " << npairs << " pairs\n";
        }

        /**********************************************************
         PHASE 2: Generate new buckets and merge clusters
        ***********************************************************/
        /* Store ids of resulting merged clusters here */
        set<int> new_clust;

        /* Populate hash table using r min-hashes */
        generate_hash_table( htable, adj, r);

        int bucketid = 1;
        /* Iterate over buckets with more than one node */
        for(auto it = htable.cbegin(); it != htable.cend(); ++it )
        {
            
            if(it->second.size()>1)
            {
                int lmin = INT_MAX -1 ;

                //cout << "\nbucket contains: ";

                for(size_t i=0; i<it->second.size(); i++)
                {
                    //cout << " " << it->second[i];

                    /* This is to distinguish nodes in same bucket later */
                    nid2bucket[it->second[i]] = bucketid;

#if(!CHECK_BUCKETS)
                    /* Add all edges between nodes in the same bucket */
                    for(size_t j=i+1; j<it->second.size(); j++)
                    {
                        gsim[it->second[i]].insert(it->second[j]); 
                        gsim[it->second[j]].insert(it->second[i]);                         
                    }          
#endif

                    /* Find minimum cid in bucket nodes */
                    lmin = ( nid2cid[it->second[i]].second < lmin ) ? nid2cid[it->second[i]].second : lmin;                        
                }

                new_clust.insert(lmin);

                bucketid++;          

                /* MERGE clusters that intersect with bucket
                    All clusters will be merged INTO the one with minimum label */
                for(auto v:it->second)
                {
                    /* Check if node cluster exists and is not C(lmin) */
                    int l_v = nid2cid[v].second;
                    //if( cid2clust.count(l_v)>0 && l_v != lmin )
                    if( l_v != lmin )
                    {
                        /* Merge nodes into C(lmin) and update their cluster labels to lmin */
                        for(auto z : cid2clust[l_v].dnode)
                        {
                            /* Update node cluster label map */
                            nid2cid[z.second].second = lmin;

                            /* Copy node into new cluster */
                            cid2clust[lmin].dnode.push_back(z);
                        }

                        /* Also merge the hidden pairs into C(lmin) */
                        for(auto h : cid2clust[l_v].hpair)
                            cid2clust[lmin].hpair.push_back(h);
#if ONE_CLIQUE
                        /* Update maxclique lower bound */
                        cid2clust[lmin].lmaxq = ( cid2clust[l_v].lmaxq > cid2clust[lmin].lmaxq ) ?
                                                  cid2clust[l_v].lmaxq : cid2clust[lmin].lmaxq ;                            
#endif                        
                        /* Delete cluster that has been merged */
                        cid2clust.erase(l_v);
                    }
                }
            /* Do NOT take bucket size as potential lmaq!
               Instead, allow it to be (possibly) detected by the max-clique.. */
            }    
        }      

        /* Clear hash table for next iteration */
        htable.clear();

        /* For every NEW cluster -> Sort nodes  wrt degrees */
        for(auto& cid: new_clust)
            if(cid2clust.count(cid)>0)
                sort(cid2clust[cid].dnode.begin(), cid2clust[cid].dnode.end(),
                     [](auto &left, auto &right) {return left.first < right.first ;});

        /**********************************************************
         PHASE 3: Add false negative edges + Fine new maximal cliques
        ***********************************************************/

        /* For every NEW cluster */
        for(auto& cid: new_clust)
        {
            if(cid2clust.count(cid)>0)
            {
                for(size_t i=0; i<cid2clust[cid].dnode.size(); i++)
                {
                    /* This is the maximum degree that satisfies min threshold */
                    size_t max_deg = (size_t) ceil( (double)cid2clust[cid].dnode[i].first /min_thres );
                    int node_a = cid2clust[cid].dnode[i].second;
                    
                    /* -> Check valid pairs to add edges / new hidden pairs */
                    size_t j = i+1;
                    while( j< cid2clust[cid].dnode.size() && cid2clust[cid].dnode[j].first <=max_deg)
                    {
                        int node_b = cid2clust[cid].dnode[j].second; 

#if(CHECK_BUCKETS)
                        bool cond = nid2cid[node_a].first != nid2cid[node_b].first; 
#else
                        bool cond = nid2cid[node_a].first != nid2cid[node_b].first &&
                                    (nid2bucket[node_a] != nid2bucket[node_b] || (nid2bucket[node_a]==0 && nid2bucket[node_b]==0) ) ; 
#endif
            
                        if(cond)
                        {
                            /* Actually compute Jaccard similarity */
                            double s = jaccard(adj, node_a, node_b);
                            //cout << "\na " << node_a<< " b " << node_b << " jac: " << s;
                            /* Put either in gsim or hidden pairs */
                            if( s >= threshold )
                            {
                                gsim[node_a].insert(node_b);
                                gsim[node_b].insert(node_a);
                            }
                            else
                            {
                                cid2clust[cid].hpair.push_back(pair<double, pair<int,int>>( s, pair<int,int>(node_a,node_b)));
                            }
                        }
                        j++;
                    }

                }
                /* Sort hidden pairs for this new cluster */ 
                sort(cid2clust[cid].hpair.begin(), cid2clust[cid].hpair.end(),
                     [](auto &left, auto &right) {return left.first < right.first;});                   
            }
        }

        /* Print clusters */
#if false
        for(auto it= cid2clust.begin(); it!= cid2clust.end(); ++it)
        {
            cout << "\nCluster: " ;
            for(auto& v: it->second.dnode)
                cout << v.second << " ";
        }
#endif

        /* For EVERY cluster look for new candidate set */
        for(auto it= cid2clust.begin(); it!= cid2clust.end(); ++it)
        {
            if(it->second.dnode.size()>1)
            {
                /* This will store canidate clique nodes */
                set<int> U;;

                /* For temporary storage of discovered cliques/ candidate sets */
                deque<int> Q;

                /* -> Check if larger clique has emerged: 
                   Iterate over every node in cluster, and check its neighborhood */
                for(size_t i=0; i<it->second.dnode.size(); i++)
                {
                    int v = it->second.dnode[i].second;
#if ONE_CLIQUE                    
                    size_t vmaxq = it->second.lmaxq;
#else
                    size_t vmaxq = node_maxq[v];
#endif
                    
                    /* Filter 1 : degree */
                    if(gsim[v].size() < vmaxq )
                        continue; 

                    /* Gather valid neighbors */
                    U.insert(v);
                    for(auto u : gsim[v])
                    {
                        if( gsim[u].size() >= vmaxq )
                        {
                            U.insert(u);
                        }
                    }

                    /* Filter 2: size of candidate clique */
                    if(U.size() <= vmaxq )
                    {
                        U.clear();
                        continue;
                    }

                    /* Check Candidate clique for maximum clique
                       This step has to be heuristic for low complexity */
                    //////////////////////////////////////////////////////////////////////
                    /* CLIQUE HEURISTIC */
                    clique_heuristic(gsim, U, Q, vmaxq);
                    if(Q.size() > vmaxq)
                    {
                        /* New max clique found! */
#if ONE_CLIQUE                        
                        it->second.lmaxq = Q.size();
#else
                        for(auto z : Q)
                            node_maxq[z] = Q.size() > node_maxq[z] ? Q.size() : node_maxq[z];
#endif

                        cout << "\n\nNew candidate: ";
                        for(auto c: Q)
                           cout << " " << c;

                        /* Add it to the candidate sets */
                        res.push_back( pair<pair<size_t, double>,set<int>>( pair<size_t, double>(Q.size(), threshold),
                                       set<int>(Q.begin(),Q.end()) ) );


                    }
                    

                    /////////////////////////////////////////////////////////////////////
                    Q.clear();
                    U.clear();
                }
            }
        }

        /* Clear bucket map */
        for(auto it = nid2bucket.begin(); it != nid2bucket.end(); ++it)
            it->second = 0;

        /* Set old_node_cluster_ids = new_node_cluster_ids  */
        for(auto it = nid2cid.begin(); it != nid2cid.end(); ++it)
            it->second.first = it->second.second;
    }
    
    /////////////////////////////////////////////////////////////////////////////////

    cout << "\nFinished generating candidates..\n";

    return res;
}




/*
Merges forward and reverse lists.
Sorts sets according to size/similarity.
Returns top-K candidate sets
*/
deque<pair<pair<size_t, double>,set<int>>> filter_candidates(   deque<pair<pair<size_t, double>,set<int>>>& cand_list,
                                                                unordered_map<int,vector<int>>& adj_fwdrev,
                                                                unordered_map<int,vector<int>>& adj_fwd,
                                                                unordered_map<int,vector<int>>& adj_rev,                                          
                                                                size_t K)
{


    //////////////////////////////////////////////////////////////////////////////////////////////////
    /* TODO: Directly detect sink/source/disc nodes from graph since they will no longer be included 
             in the corresponding hash-tables*/
    //////////////////////////////////////////////////////////////////////////////////////////////////
    set<int> sinks;
    set<int> sources;
    set<int> disc;

    /* Put aside sink nodes (do not point anywhere) */
    for(auto it=adj_fwd.cbegin(); it != adj_fwd.cend(); ++it)
        if(it->second.empty())
            sinks.insert(it->first);
            
     /* Put aside source nodes (are not pointed) */       
    for(auto it=adj_rev.cbegin(); it != adj_rev.cend(); ++it)
        if(it->second.empty())
            sources.insert(it->first);

    /* Find completely incative nodes (neither point, nor pointed, SELF-LOOPs possible)  */        
    set_intersection(sources.begin(), sources.end(), sinks.begin(), sinks.end(), inserter( disc, disc.begin()));

    /* Clean sinks and sources from completely disconected nodes */
    for(auto& v : disc)
    {
        sinks.erase(v);
        sources.erase(v);
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////

    /* Check for potential INSTAR spokes
       This will be nodes with similar outgoing connections and little-to-no incoming connections
       If they have 1 common destination, add it to the set (since it will be the hub)*/
    set<int> temp;   
    for(auto it = cand_list.begin(); it != cand_list.end(); ++it)
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
                    if(hub != SELF_ID)
                    {
                        it->second.insert(hub);
                        it->first.first++;
                    }
                    // cout << "\n\n Found INSTAR with hub: " << hub
                    //      << ", size: " << it->second.size() << ", score " << it->first.second << "\n\n"; 
                } 
            }

            //cout << "\nrevlist " << it->second.size() << " sources " << sources.size() 

            set_intersection(it->second.begin(), it->second.end(), sinks.begin(),
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
                        if(hub != SELF_ID)
                        {
                            it->second.insert(hub);
                            it->first.first++;
                            // cout << "\n\n Found OUTSTAR with hub: " << hub
                            //      << ", size: " << it->second.size() << ", score " << it->first.second << "\n\n";                         

                        }
                    } 
                }
            }
        }              
           
    }

        
    /* Merge lists */
    deque<pair<pair<size_t, double>,set<int>>> res;
    for(auto f : cand_list)
        res.push_back( f );       

#if(BUCKETIZE_SELF_LOOPS)
    set<int> sl;
    for(auto it=adj_fwd.cbegin(); it != adj_fwd.cend(); ++it)
        if(it->second.size()==1)
            if( it->first == it->second[0])
                sl.insert(it->first);
    res.push_back( pair<pair<size_t, double>,set<int>>( pair<size_t, double>(sl.size(),1.0), sl ));    
#endif

#if(BUCKETIZE_EXCHANGE_PAIRS)
    for(auto it=adj_fwd.cbegin(); it != adj_fwd.cend(); ++it)
        if(it->second.size()==1)
            if(adj_fwd[it->second[0]].size()==1)
                if( it->first == adj_fwd[it->second[0]][0] && it->first != it->second[0] )
                {
                    set<int> ex_pair {it->first,it->second[0]};
                    res.push_back( pair<pair<size_t, double>,set<int>>( pair<size_t, double>(2,1.0), ex_pair)); 
                }    
#endif


    /* Sort */
    sort(res.begin(), res.end(), 
         [](auto &left, auto &right) {return ((double) left.first.first)*left.first.second > ((double) right.first.first)*right.first.second;} );
    // sort(res.begin(), res.end(), 
    //      [](auto &left, auto &right) {return ((double) left.first.first) > ((double) right.first.first);} );
    /* Add disconnected set in the beginning since it will almost definetly help reduce cost */
    res.push_front( pair<pair<size_t, double>,set<int>>( pair<size_t, double>(disc.size(),1.0), set<int>(disc.begin(),disc.end())));

    /* Return only the best K  */
    while(res.size() > K)
        res.pop_back();
    return res;
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////

/*------------------------------------------------------------------------

                        SUMMARIZATION ALGORITHM(S)

-------------------------------------------------------------------------*/


/////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
 Measure cost reduction from candidate set merging 
*/
double cost_reduction( double cost_before, graph& graph_,  set<int>& cand_set )
{
    if(cand_set.size() <= 1)
        return 0;

#if PROFILE_MRG_EXP_SUB
    struct timeval start = {.tv_sec = 0, .tv_usec = 0 };
    struct timeval stop = {.tv_sec = 0, .tv_usec = 0 };
#endif

    /* Separate candidate set if types dont agree */
    unordered_map< int, vector<int> > cand_per_type;
    for(auto& v : cand_set)
    {
        if(cand_per_type.count(graph_.snode[v].type) == 0 )
            cand_per_type.insert( pair<int, vector<int>>( graph_.snode[v].type, vector<int>()));
        cand_per_type[graph_.snode[v].type].reserve(cand_per_type[graph_.snode[v].type].size() + 1 );
        cand_per_type[graph_.snode[v].type].push_back(v);
    }

#if PROFILE_MRG_EXP_SUB
    gettimeofday( &start, NULL);
#endif

    /* Merge sets and save ids of new super-nodes */
    vector<int> new_supernodes;
    new_supernodes.reserve(cand_per_type.size());
    for(auto it = cand_per_type.begin(); it != cand_per_type.end(); ++it)
        if(it->second.size() > 1)
            new_supernodes.push_back(graph_.merge_snodes( it->second ));

    /* Compute cost after merging */
    pair<double,double> costs_after = graph_.get_cost();
    double cost_after = costs_after.first + costs_after.second;

#if PROFILE_MRG_EXP_SUB
    gettimeofday( &stop, NULL);
    auto dur =  stop.tv_usec - start.tv_usec;
    cout << "\n\nmrg dur " << (size_t) dur;
#endif

    //graph_.cleanup();

#if PROFILE_MRG_EXP_SUB
    gettimeofday( &start, NULL);
#endif

    /* Undo merging by expanding new super-nodes */
    graph_.expand_snodes(new_supernodes);    

#if PROFILE_MRG_EXP_SUB
    gettimeofday( &stop, NULL);
    dur =  stop.tv_usec - start.tv_usec;
    cout << " exp dur " << (size_t) dur;

    gettimeofday( &start, NULL);
#endif

    graph_.cleanup();

    pair<double,double> c = graph_.get_cost();
    double ct = c.first + c.second;

#if PROFILE_MRG_EXP_SUB
    gettimeofday( &stop, NULL);
    dur =  stop.tv_usec - start.tv_usec;
    cout << " cleanup_cost dur " << (size_t) dur;
#endif


    //cout << "\nbefore " << cost_before << " rec " << ct << "\n";
    assert( fabs(ct - cost_before) <= EPS * fabs(ct) );

    return cost_before - cost_after;
}



/*
    Works by finding weighted independent set over the graph of overlapping candidates
    The weight of every node (candidate) is equal to its ``delta'' value (computed once)
    Returns total cost before and after summarization.
*/
pair<double, double> wis_based_summ(graph& graph_, unordered_map<int,vector<int>>& adj,
                                    deque<pair<pair<size_t, double>,set<int>>>& candidates, string save_all_to)
{
    pair<double,double> res;

    unordered_map< int, vector<int> > cand_per_type;

    size_t iter_ = 0;

    unordered_map<int, deque<int>> inverted_ind = cand_set_invert_index(candidates);
    unordered_map<int, unordered_set<int>> overlap_graph = get_cand_overlap_graph(inverted_ind, candidates.size());
    

    vector<double> delta( candidates.size(), 0.0);

    if(save_all_to.compare("IGNORE"))
        graph_.save_as(save_all_to + "/" + to_string(iter_));

    pair<double,double> costs_before = graph_.get_cost();
    double cost_before = costs_before.first + costs_before.second;
    res.first = cost_before ;

#if true
    /* Do one pass over all candidates to compute all deltas */
    cout << "\n\nFirst pass.. Computing deltas..\n\n";
    for(size_t i=0; i< candidates.size(); i++)
    {
        delta[i] = cost_reduction( cost_before, graph_, candidates[i].second );
        cout <<"\n\n"<< i+1 <<" out of " << candidates.size(); 
        cout << "\n delta = " << delta[i];
    }
#else
    for(size_t i=0; i< candidates.size(); i++)
        delta[i] = candidates[i].first.second * (double) candidates[i].first.first;
#endif

    cout << "\n\nFinding WIS candidate sets..\n";

    /* Use deltas as node weights and approximate the weighted indepenent
       set over the candidate-set-overlap graph  */
    double delta_sum = 0.0;   
    unordered_set<size_t> best_cands = approx_wis( delta, overlap_graph, delta_sum);    

    cout << "\n\nMerging " << best_cands.size() << " WIS candidate sets..\n";

    for(auto best_cand: best_cands)
    {
        /* MERGE BEST CANDIDATE SET */
        ///////////////////////////////////////////////////////////
        // cout << "\nBest Candidate: ";
        // for(auto& v : candidates[best_cand].second) cout << v << " ";
        // cout << "\n";

        /* Separate candidate set if types dont agree */
        for(auto& v : candidates[best_cand].second)
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
    }

    /* Get cost of graph after summarization */
    pair<double,double> costs_after = graph_.get_cost();
    double cost_after = costs_after.first + costs_after.second;

    res.second = cost_after;

    cout << "\n\nDelta_sum " << delta_sum << ", True delta " << res.first - res.second << "\n"; 

    return res;

}


/*
    Greedy summarization of graph given list of candidate sets.
    Summarizes graph in place.
    Returns total cost before and after summarization.
*/
pair<double, double> more_greedy_light_summ(graph& graph_, unordered_map<int,vector<int>>& adj,
                                            deque<pair<pair<size_t, double>,set<int>>>& candidates, string save_all_to)
{
    pair<double,double> res;

    unordered_map< int, vector<int> > cand_per_type;

    size_t iter_ = 0;


    unordered_map<int, deque<int>> inverted_ind = cand_set_invert_index(candidates);
    unordered_map<int, unordered_set<int>> overlap_graph = get_cand_overlap_graph(inverted_ind, candidates.size());
    unordered_map<int, unordered_set<int>> prox_graph = get_cand_prox_graph(inverted_ind, adj, candidates.size());  
    

    vector<bool> valid_cand( candidates.size(), true);
    vector<double> delta( candidates.size(), 0.0);
    set<int> to_be_updated;

    if(save_all_to.compare("IGNORE"))
        graph_.save_as(save_all_to + "/" + to_string(iter_));


    /* Do one pass over all candidates to compute all deltas */
    cout << "\n\nFirst pass.. Computing deltas..\n\n";
    pair<double,double> costs_before = graph_.get_cost();
    double cost_before = costs_before.first + costs_before.second;
    res.first = cost_before ;
    for(size_t i=0; i< candidates.size(); i++)
    {
        delta[i] = cost_reduction( cost_before, graph_, candidates[i].second );
        cout <<"\n"<< i+1 <<" out of " << candidates.size(); 
    }
        
 

    while(true)
    {
        cout << "\n\nIteration: " << ++iter_;

        /* Find best canidate at this iteration */
        size_t best_cand=0;
        double best_delta = 0.0;
        for(size_t i=0; i<candidates.size(); i++)
        {
            if( valid_cand[i] && delta[i] > best_delta)
            {
                best_delta = delta[i];
                best_cand = i;
            }
        }

        cout << "\n\nBest delta: " << best_delta;

        /* Terminate if no improvement possible */
        if(best_delta <= 0.0)
        {
            cout << "\n\nNo further improvement possible..\n\n";
            break; 
        }

        /* MERGE BEST CANDIDATE SET */
        ///////////////////////////////////////////////////////////
        cout << "\nBest Candidate: ";
        for(auto& v : candidates[best_cand].second) cout << v << " ";
        cout << "\n";
        /* Separate candidate set if types dont agree */
        for(auto& v : candidates[best_cand].second)
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

        /* Save summary graph */
        if(save_all_to.compare("IGNORE"))
            graph_.save_as(save_all_to + "/" + to_string(iter_));
        /////////////////////////////////////////////////////////////

        /* UPDATE CANDIDATE GRAPHS AND DELTAS */


        /* Update ovelapping sets */
        valid_cand[best_cand] = false;
        unordered_set<int> this_overlap = overlap_graph[best_cand];
        for(auto i : this_overlap)
        {

            if(valid_cand[i] )
            {
                /* Faster set difference in place */
                auto a_begin = candidates[i].second.begin();
                auto a_end = candidates[i].second.end();
                auto b_begin = candidates[best_cand].second.begin();
                auto b_end = candidates[best_cand].second.end();
                 
                for (auto ita = a_begin, itb = b_begin; ita != a_end && itb != b_end;) {
                    if (*ita < *itb) {
                        ++ita;
                    } else if (*itb < *ita) {
                        ++itb;
                    } else {
                        ita = candidates[i].second.erase(ita);
                        ++itb;
                    }
                }
            }
            
            if(candidates[i].second.size() <= 1)
            {
                valid_cand[best_cand] = false;
                cand_graph_remove_node(overlap_graph, i);
                cand_graph_remove_node(prox_graph, i);
            }

            to_be_updated.insert(i);
        }
        cand_graph_remove_node(overlap_graph, best_cand);



        /* Check neighboing sets */
        for(auto i : prox_graph[best_cand])
            to_be_updated.insert(i);
        cand_graph_remove_node(prox_graph, best_cand);

        /* Update deltas of neighbors */
        costs_before = graph_.get_cost();
        cost_before = costs_before.first + costs_before.second;
        size_t updt_cnt = 0;
        for(auto i : to_be_updated)
            if(valid_cand[i])
            {
                delta[(size_t) i] = cost_reduction( cost_before, graph_, candidates[ (size_t) i].second );
                updt_cnt++;
                cout << "\n\nupd " << updt_cnt << " size " << candidates[ (size_t) i].second.size();
            }
        cout << "\n\nDeltas updated: " << updt_cnt << "\n";        

        to_be_updated.clear();

        /* Save summary graph */
        if(save_all_to.compare("IGNORE"))
            graph_.save_as(save_all_to + "/" + to_string(iter_));

    }

    /* Get cost of graph after summarization */
    pair<double,double> costs_after = graph_.get_cost();
    double cost_after = costs_after.first + costs_after.second;

    res.second = cost_after;

    return res;

}



/*
    Greedy summarization of graph given list of candidate sets.
    Summarizes graph in place.
    Returns total cost before and after summarization.
*/
pair<double, double> more_greedy_summ(graph& graph_,  deque<pair<pair<size_t, double>,set<int>>>& candidates, 
                                      string save_all_to)
{
    pair<double,double> res;

    unordered_map< int, vector<int> > cand_per_type;

    size_t iter_ = 0;

    vector<bool> valid_cand( candidates.size(), true);

    if(save_all_to.compare("IGNORE"))
        graph_.save_as(save_all_to + "/" + to_string(iter_));

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
            if(valid_cand[cand_counter])
            {
                double delta = cost_reduction( cost_before, graph_, cand_set.second );
                cout << "\n delta: " <<  delta <<"\n";
                if(delta > best_delta)
                {
                    best_delta = delta;
                    best_cand = cand_counter;
                }
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
        for(auto& v : candidates[best_cand].second) cout << v << " ";
        cout << "\n";
        /* Separate candidate set if types dont agree */
        for(auto& v : candidates[best_cand].second)
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

        /* Save summary graph */
        if(save_all_to.compare("IGNORE"))
            graph_.save_as(save_all_to + "/" + to_string(iter_));

        /* NOTE: Bellow operations may be slow!*/

        /* Subtract nodes that were summarized from other candidate sets */
        for(size_t i=0; i<candidates.size(); i++)
        {
            if(i!= best_cand && valid_cand[i] )
            {
                /* Faster set difference in place */
                auto a_begin = candidates[i].second.begin();
                auto a_end = candidates[i].second.end();
                auto b_begin = candidates[best_cand].second.begin();
                auto b_end = candidates[best_cand].second.end();
                 
                for (auto ita = a_begin, itb = b_begin; ita != a_end && itb != b_end;) {
                    if (*ita < *itb) {
                        ++ita;
                    } else if (*itb < *ita) {
                        ++itb;
                    } else {
                        ita = candidates[i].second.erase(ita);
                        ++itb;
                    }
                }
            }
        }


        /* Remove the candidate set that was used, and empty/singleton sets */
        valid_cand[best_cand] = false;
        for(size_t i=0; i<candidates.size(); i++)
            if( candidates[i].second.size() <= 1)
                valid_cand[i] = false;
                

        /* Also terminate if set becomes empty */
        if(!any_of(valid_cand.begin(), valid_cand.end(), [](bool x) {return x;} ))
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



/*
    Greedy summarization of graph given list of candidate sets.
    Summarizes graph in place.
    Returns total cost before and after summarization.
*/

pair<double, double> greedy_summ(graph& graph_,  deque<pair<pair<size_t, double>,
                                 set<int>>>& candidates, string save_all_to)
{
    pair<double,double> res;

    unordered_map< int, vector<int> > cand_per_type;

    size_t iter_ = 0;

    set<int> merged_set, merged_snodes;

    /* Save summary graph */
    if(save_all_to.compare("IGNORE"))
        graph_.save_as(save_all_to + "/" + to_string(iter_));    

    while(iter_ < candidates.size())
    {
        cout << "\n\nIteration: " << +iter_;

        /* Subtract nodes that have been summarized candidate set */
        auto a_begin = candidates[iter_].second.begin();
        auto a_end = candidates[iter_].second.end();
        auto b_begin = merged_set.begin();
        auto b_end = merged_set.end();
            
        for (auto ita = a_begin, itb = b_begin; ita != a_end && itb != b_end;) {
            if (*ita < *itb) {
                ++ita;
            } else if (*itb < *ita) {
                ++itb;
            } else {
                ita = candidates[iter_].second.erase(ita);
                ++itb;
            }
        }


        /* Get cost before merging */
        pair<double,double> costs_before = graph_.get_cost();
        double cost_before = costs_before.first + costs_before.second;

        if(iter_ == 1)
            res.first = cost_before;

        /* MERGE CANDIDATE */
        ///////////////////////////////////////////////////////////

        /* Separate candidate set if types dont agree */
        for(auto& v : candidates[iter_].second)
        {
            if(cand_per_type.count(graph_.snode[v].type) == 0 )
                cand_per_type.insert( pair<int, vector<int>>( graph_.snode[v].type, vector<int>()));
            cand_per_type[graph_.snode[v].type].reserve(cand_per_type[graph_.snode[v].type].size() + 1 );
            cand_per_type[graph_.snode[v].type].push_back(v);
        }
            
        /* Merge sets per type */
        for(auto it = cand_per_type.begin(); it != cand_per_type.end(); ++it)
            if(it->second.size() > 1)
                merged_snodes.insert(graph_.merge_snodes( it->second ));
            

        graph_.cleanup();
        cand_per_type.clear();        
        /////////////////////////////////////////////////////////////

        /* Get cost after merging */
        pair<double,double> costs_after = graph_.get_cost();
        double cost_after = costs_after.first + costs_after.second;                

        /* If cost did not reduce then expand again */
        double delta = cost_before - cost_after;
        cout << "\n\n    Size: " << candidates[iter_].second.size();
        cout << "\n     Sim: " << candidates[iter_].first.second;
        cout << "\n     Delta: " << delta;
        if(delta <= 0.0)
        {
                vector<int> v = vector<int>(merged_snodes.cbegin(), merged_snodes.cend());
                graph_.expand_snodes(v);
        }
        else
        {
                /* Keep the merging */
                merged_set.insert(candidates[iter_].second.begin(), candidates[iter_].second.end());
        
                /* Save summary graph */
                if(save_all_to.compare("IGNORE"))
                    graph_.save_as(save_all_to + "/" + to_string(iter_+1));
        }
        
        //graph_.cleanup();

        merged_snodes.clear();
        iter_++;
    }

    /* Get cost of graph after summarization */
    pair<double,double> costs_after = graph_.get_cost();
    double cost_after = costs_after.first + costs_after.second;


    res.second = cost_after;

    return res;

}