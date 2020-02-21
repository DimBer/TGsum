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
std::deque<std::pair<std::pair<size_t, double>,std::set<int>>> cand_gen(std::unordered_map<int,std::vector<int>>& adj,
                                                                              size_t r, size_t b);


/*
Merges forward and reverse lists.
Sorts sets according to size/similarity.
Returns top-K candidate sets
*/
std::deque<std::pair<std::pair<size_t, double>,std::set<int>>> filter_candidates(   std::deque<std::pair<std::pair<size_t, double>,std::set<int>>>& fwd_list,
                                                                                    std::unordered_map<int,std::vector<int>>& adj_fwdrev,
                                                                                    std::unordered_map<int,std::vector<int>>& adj_fwd,
                                                                                    std::unordered_map<int,std::vector<int>>& adj_rev,
                                                                                    size_t K);


/* 
Measure cost reduction from candidate set merging 
*/
double cost_reduction( double cost_before, graph& graph_,  std::set<int>& cand_set );


/*
    Very Greedy summarization of graph given list of candidate sets.
    Summarizes graph in place.
    Returns total cost before and after summarization.
*/
std::pair<double, double> more_greedy_summ(graph& graph_,  std::deque<std::pair<std::pair<size_t, 
                                           double>,std::set<int>>>& candidates, std::string save_all_to);

/*
    Very Greedy summarization of graph given list of candidate sets.
    Summarizes graph in place.
    Graph assisted to save operations in large
    Returns total cost before and after summarization.
*/
std::pair<double, double> more_greedy_light_summ(graph& graph_,  std::unordered_map<int,std::vector<int>>& adj,
                                           std::deque<std::pair<std::pair<size_t, double>,std::set<int>>>& candidates, std::string save_all_to);

/*
    Works by finding weighted independent set over the graph of overlapping candidates
    The weight of every node (candidate) is equal to its ``delta'' value (computed once)
    Returns total cost before and after summarization.
*/
std::pair<double, double> wis_based_summ(graph& graph_,  std::unordered_map<int,std::vector<int>>& adj,
                                         std::deque<std::pair<std::pair<size_t, double>,std::set<int>>>& candidates, std::string save_all_to);


/*
    Very Greedy summarization of graph given list of candidate sets.
    Summarizes graph in place.
    Returns total cost before and after summarization.
*/
std::pair<double, double> greedy_summ(graph& graph_,  std::deque<std::pair<std::pair<size_t,
                                      double>,std::set<int>>>& candidates, std::string save_all_to);


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