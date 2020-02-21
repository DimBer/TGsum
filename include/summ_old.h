#ifndef SUMM_H
#define SUMM_H

#include <vector>
#include <memory>
#include <boost/functional/hash.hpp>

#include "graph.h"

/*
    Hash (LSH) nodes into sets (sorted vectors).
    Return ranked list of sets.  
*/
std::vector<std::pair<std::pair<size_t, double>,std::set<int>>> increment_lsh(std::unordered_map<int,std::vector<int>>& adj,
                                                                              size_t r, size_t b);


/*
Merges forward and reverse lists.
Sorts sets according to size/similarity.
Returns top-K candidate sets
*/
std::vector<std::set<int>> get_candidates(std::vector<std::pair<std::pair<size_t, double>,std::set<int>>>& fwd_list,
                                          std::vector<std::pair<std::pair<size_t, double>,std::set<int>>>& rev_list,
                                          std::unordered_map<int,std::vector<int>>& adj_fwd,
                                          std::unordered_map<int,std::vector<int>>& adj_rev,
                                          size_t K);


/* 
Measure cost reduction from candidate set merging 
*/
double cost_reduction( double cost_before, graph& graph_,  std::set<int>& cand_set );


/*
    Greedy summarization of graph given list of candidate sets.
    Summarizes graph in place.
    Returns total cost before and after summarization.
*/
std::pair<double, double> greedy_summ(graph& graph_,  std::vector<std::set<int>>& candidates);


/*
    Template wrapper for BOOST generic range_based hasher 
*/
template <typename Container> 
struct container_hash {
    std::size_t operator()(Container const& c) const {
        return boost::hash_range(c.begin(), c.end());
    }
};


#endif