#ifndef UTILS_H
#define UTILS_H

#define SELF_ID 01020304

#include <vector>
#include <set>
#include <unordered_set>
#include <deque>
#include <unordered_map>

/*
Types of glyphs
*/
enum glyph_type {DISC, CLIQUE, INSTAR, OUTSTAR, NONE};

/*
Dichotomous search over multiplicity encoding cost
*/
std::pair<size_t, double> dichotomous_search(const std::vector<size_t>& mults);


/*
Binomial encoding cost os subset of k from n: L_N(k) + log_2(n choose k)
*/
double binom_cost(size_t n, size_t k);

/*
Returns ration of false positive LSH pairs over all possible pairs 
*/
double check_FP( const std::set<int>& nodes, std::unordered_map<int, std::vector<int>>& adj, double threshold);

/*
 Returns true average similarity of nodes in cluster
*/
double  get_sim(const std::set<int>& nodes, std::unordered_map<int, std::vector<int>>& adj); 


/* Intersection and union size of SORTED vectors */
size_t union_size(std::vector<int>& vec_a, std::vector<int>& vec_b);

size_t intersection_size(std::vector<int>& vec_a, std::vector<int>& vec_b);

/* Save candidate sets to file */
void store_candidates(std::string cands_file_path, std::deque<std::pair<std::pair<size_t, double>,std::set<int>>>& candidates);


/* Delete all files in directory */
void empty_dir(std::string);

/*
  Save the adjacency lists to file
*/
void save_cand_graph(std::unordered_map<int, std::unordered_set<int>>& adj, std::string filename);


/*
  ''Splice'' the forawrd and backward adjacency lists into a single forward/backward representation
    To preserve the set property, when v_to = {a ,b ,c}, and v_from = {d, e, f}, the total representation will be
    v_tofrom = {a,b,c, d + v_max, e + v_max, f + v_max}.
*/
std::unordered_map<int,std::vector<int>> get_forward_backward_lists(std::unordered_map<int,std::vector<int>>& fwd_graph,
                                                                   std::unordered_map<int,std::vector<int>>& rev_graph);


/*
  Ratio of nodes covered by candidate sets
*/
double get_cand_coverage(size_t N, std::deque<std::pair<std::pair<size_t, double>,std::set<int>>>& candidates);


#endif