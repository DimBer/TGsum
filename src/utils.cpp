#include <cstdlib>
#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <random>
#include <cassert>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <deque>
#include <dirent.h>
#include "omp.h"

#include "utils.h"


using namespace std;

/*------------------------------------------------------------------------

                        Frequently used functions

-------------------------------------------------------------------------*/

/*
Binomial encoding cost of encoding subset of k from n: L_N(k) + log_2(n choose k)
*/
double binom_cost(size_t n, size_t k)
{
    double res = 0.0;

    if(n == 0 && k == 0)
        return 0.0;

    if(k == 0)
        return 0.0;

    if(k == n)
        return 2.0* log2( (double) k) + 1.0;

    res += 2.0* log2( (double) k) + 1.0;

    for(size_t i = n-k+1; i <= n; i++)
        res += log2( (double) i);

    for(size_t i = 1; i <= k; i++)
        res -= log2( (double) i);

    return res;
}




/*
Encoding cost
*/
static double get_cost(const vector<size_t>& vm, double x_)
{
    double res = 0.0;
    size_t x = (size_t) x_;

    for(auto& m : vm)
        res += (x == m) ? 1.0 : 2.0* log2(abs((double) ( (int)m-(int)x))) + 2.0 ;

    return res;
}

/*
Dichotomous search over multiplicity encoding cost
*/
pair<size_t, double> dichotomous_search(const vector<size_t>& mults)
{

//    cout << "MULTS: " ; 
//    for(auto x: mults) 
//       cout << " " << x; 
//    cout << "\n" ; 

   double x_r = *max_element(mults.begin(), mults.end());
   double x_l = *min_element(mults.begin(), mults.end());
   double x_m = x_l;
   double cost_l = get_cost(mults, x_l);
   double cost_r = get_cost(mults,x_r);


   while( x_l != x_r )
   {
       x_m = (x_l + x_r)/2.0;
       if(cost_l < cost_r)
       {       
           x_r = x_m;
           cost_r = get_cost(mults, x_r);
       }
       else
       {
           x_l = x_m;
           cost_l = get_cost(mults, x_l);
       }
    //    cout << " left " << x_l << " right " << x_r <<"\n";
    //    cout << " left " << cost_l << " right " << cost_r <<"\n";
   };

//    cout << "Rep: " << x_m; 

   return pair<size_t, double>((size_t) x_l , cost_l);

}




/*
Returns ratio of false positive LSH pairs over all possible pairs 
*/
double check_FP( const set<int>& nodes_, unordered_map<int,vector<int>>& adj, double threshold)
{

    vector<int> nodes(nodes_.cbegin(), nodes_.cend());

    double all_pairs = nodes.size()*(nodes.size() - 1)/2; 
    size_t fp_pairs = 0;

    vector<int> temp;
    temp.reserve(adj.size());

    for(size_t i=0; i<nodes.size(); i++)
    {
        for(size_t j=i+1; j<nodes.size(); j++)
        {
            auto it3 = set_intersection(  adj[nodes[i]].cbegin(), adj[nodes[i]].cend(), adj[nodes[j]].cbegin(), adj[nodes[j]].cend(), temp.begin());
            temp.resize(it3-temp.begin()); 
            double intersection = (double) temp.size();

            auto it4 = set_union(  adj[nodes[i]].cbegin(), adj[nodes[i]].cend(), adj[nodes[j]].cbegin(), adj[nodes[j]].cend(), temp.begin());
            temp.resize(it4-temp.begin()); 
            double union_ = (double) temp.size();

            if( intersection / union_ < threshold )
                fp_pairs ++;

        }
    }


    return (double) fp_pairs/ all_pairs;
}





/*
 Returns true average similarity of nodes in cluster
*/
double get_sim( const set<int>& nodes_, unordered_map<int,vector<int>>& adj)
{

    vector<int> nodes(nodes_.cbegin(), nodes_.cend());

    double all_pairs = nodes.size()*(nodes.size() - 1)/2; 
    double res = 0.0;

    for(size_t i=0; i<nodes.size(); i++)
    {
        for(size_t j=i+1; j<nodes.size(); j++)
        {
            double intersection = (double) intersection_size(adj[nodes[i]], adj[nodes[j]]);

            double union_ = (double) union_size(adj[nodes[i]], adj[nodes[j]]);

            res += intersection/union_;
        }
    }


    return res/ all_pairs;
}


/* Measure intersection size of two SORTED vectors */
size_t intersection_size(vector<int>& vec_a, vector<int>& vec_b)
{
    if(vec_a.size() == 0 || vec_b.size() == 0)
        return 0;

    size_t count = 0, i =0, j = 0;

    while( i < vec_a.size() && j < vec_b.size() )
    {
        if(vec_a[i] == vec_b[j])
        {
            count++;
            i++;
            j++;
        }
        else
        {
            if(vec_a[i] > vec_b[j])
            {
                j++;
            }
            else
            {
                i++;
            }
        }
    }

    return count;
}


/* Measure union size of two sorted ivec's */
size_t union_size(vector<int>& vec_a, vector<int>& vec_b)
{
    if(vec_a.size() == 0 && vec_b.size() == 0)
        return 0;

    if(vec_a.size() == 0 && vec_b.size() > 0)
        return vec_b.size();

    if(vec_a.size() > 0 && vec_b.size() == 0)
        return vec_a.size();

    size_t count = 0, i =0, j = 0;

    while( i < vec_a.size() && j < vec_b.size() )
    {
        if(vec_a[i] == vec_b[j])
        {
            i++;
            j++;
        }
        else
        {
            if(vec_a[i] > vec_b[j])
            {
                j++;
            }
            else
            {
                i++;
            }
        }
        count++;
    }

    while(i<vec_a.size())
    {
        i++;
        count++;
    } 

    while(j<vec_b.size())
    {
        j++;
        count++;
    }

    return count;
}




/* Save candidate sets to file */
void store_candidates(string cands_file_path, deque<pair<pair<size_t, double>, set<int>>>& candidates)
{

    ofstream myfile;
    myfile.open(cands_file_path);

    for(auto cand : candidates)
    {
        myfile << cand.first.second << " ";
        for(auto v : cand.second)
            myfile << v << " ";
        myfile << "\n";
    }

    myfile.close();
}


/* Delete all files in directory */
void empty_dir(string dir_path)
{
  
    string dir_path_ = dir_path + "/";
    const char* cdir_path = dir_path_.c_str();

    DIR *theFolder = opendir(cdir_path);
    struct dirent *next_file;
    char filepath[256];

    while ( (next_file = readdir(theFolder)) != NULL )
    {
        // build the path for each file in the folder
        sprintf(filepath, "%s/%s", cdir_path, next_file->d_name);
        remove(filepath);
    }
    closedir(theFolder);
}

/*
  Save the adjacency lists to file
*/
void save_cand_graph(unordered_map<int, unordered_set<int>>& adj,string filename)
{

    ofstream myfile;
    myfile.open(filename);

    for(auto it = adj.cbegin(); it != adj.cend(); ++it)
    {
        int node_a = it->first;
        for(auto node_b : it->second)
            myfile << node_a << " " << node_b <<"\n";
    }

    myfile.close();
}


/*
  ''Splice'' the forawrd and backward adjacency lists into a single forward/backward representation
    To preserve the set property, when v_to = {a ,b ,c}, and v_from = {d, e, f}, the total representation will be
    v_tofrom = {a,b,c, d + v_max, e + v_max, f + v_max}.
*/
unordered_map<int, vector<int>> get_forward_backward_lists(unordered_map<int, vector<int>>& fwd_graph,
                                                           unordered_map<int, vector<int>>& rev_graph)
{

    unordered_map<int, vector<int>> res;

    int max_node = 0;
    for(auto it = fwd_graph.cbegin(); it != fwd_graph.cend(); ++it)
    {
        max_node = it->first > max_node ? it->first : max_node;
        res.insert( pair<int, vector<int>>( it->first, vector<int>(it->second.cbegin(), it->second.cend())));
    }

    for(auto it = rev_graph.cbegin(); it != rev_graph.cend(); ++it)
        max_node = it->first > max_node ? it->first : max_node;

    for(auto it = rev_graph.cbegin(); it != rev_graph.cend(); ++it)
    {

        vector<int> temp(it->second.size());

        for(size_t i = 0; i< it->second.size(); i++)
            temp[i] = it->second[i] + max_node;

        if(res.count(it->first)==0)
        {
            res.insert( pair<int, vector<int>>( it->first, temp) );
        }
        else
        {
            res[it->first].reserve(res[it->first].size() + temp.size() );
            for(auto v : temp)
                res[it->first].push_back(v);
        }
    }

    // for(auto it= res.begin(); it != res.end(); ++it)
    // {
    //     cout << "\nNode: " << it->first << " -> ";
    //     sort(it->second.begin(), it->second.end());
    //     for(auto& v:it->second) cout << " " << v;
    // }
        

    return res;
}

/*
  Ratio of nodes covered by candidate sets
*/
double get_cand_coverage(size_t N, deque<pair<pair<size_t, double>,set<int>>>& candidates)
{
    unordered_set<int> all;

    for(auto& cand : candidates)
        for(auto v : cand.second)
            all.insert(v);

    return (double) all.size()/ (double) N;
}