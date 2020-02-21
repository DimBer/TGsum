#include <cstdlib>
#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <random>
#include <cassert>
#include <chrono>
#include <set>
/*
DRIVER PROGRAMM TO TEST MERGING/ SUMMARIZATION
*/

#include "graph.h"
#include "summ_opts.h"
#include "summ.h"


using namespace std;


int main(int argc, char **argv)
{

  /* Get options from cmd args */
  struct my_opt opts = get_options(argc, argv);

  /* If saving intermediate results for the same graph, first delete previous ones*/
  cout << opts.save_all_to.compare("IGNORE")<<"\n"; 
  if(opts.save_all_to.compare("IGNORE"))
    empty_dir(opts.save_all_to);  

 /////////////////////////////////////////////////////////////////////////////////////////////
                    /*  TEST SUMMARIZE  */
 /////////////////////////////////////////////////////////////////////////////////////////////

  size_t nseeds = opts.niters;
  double ratio = 0.0;
  double cand_coverage = 0.0;
  size_t runtime = 0;
  graph test_graph;
  for(int seed = 1;  seed< (int) nseeds+1; seed++)
  {
    srand(seed);

    test_graph = graph(opts.graph, opts.types, opts.no_mult, opts.no_type);

    auto fwd_graph = test_graph.get_fwd_graph();
    auto rev_graph = test_graph.get_rev_graph();
    auto fwdrev_graph = get_forward_backward_lists(fwd_graph,rev_graph);

    auto start = chrono::high_resolution_clock::now(); 

    cout << "\n\nHashind forward/reverse graph..\n\n";
    auto cand_raw = cand_gen(fwdrev_graph, opts.rows, opts.bands );
    cout << "\n\nMerging lists and generating candidates\n\n";  
    auto candidates = filter_candidates(cand_raw, fwdrev_graph, fwd_graph, rev_graph, opts.K);

    string cands_file_path = "data/cand/out/r" + to_string(opts.rows) ;
    if(opts.save_cands)
      store_candidates(cands_file_path, candidates);

    cand_raw.clear();

    cout << "\n\nBegin summarization\n\n";
    pair<double,double> bna_costs;
    if(opts.wis_based)
      bna_costs = wis_based_summ(test_graph, fwd_graph, candidates, opts.save_all_to);
    else if(opts.more_greedy_light)
      bna_costs = more_greedy_light_summ(test_graph, fwd_graph, candidates, opts.save_all_to);
    else if(opts.more_greedy)
      bna_costs = more_greedy_summ(test_graph, candidates, opts.save_all_to);
    else
      bna_costs = greedy_summ( test_graph, candidates, opts.save_all_to);

    test_graph.display();

    auto stop = chrono::high_resolution_clock::now(); 
  
    auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);

    cand_coverage += get_cand_coverage(test_graph.nnodes_full, candidates);

    ratio += bna_costs.second/bna_costs.first;
    runtime += duration.count();

  }

  /* Save resulting summary graph */
  if(opts.save_to.compare("IGNORE"))
    test_graph.save_as( opts.save_to + "_final" );


  cout << "\nMean compression ratio: " << ratio/ (double) nseeds;
  cout << "\nMean runtime: " << runtime / nseeds;

  cout << "\nMean candidate coverage: " << cand_coverage/ (double) nseeds;

  cout << "\n\nDONE\n\n";







 /////////////////////////////////////////////////////////////////////////////////////////////
                    /*  TEST EXPANDING  */
 ////////////////////////////////////////////////////////////////////////////////////////////

  // pair<double,double> cost;

  // graph test_graph = graph(opts.graph, opts.types, opts.no_mult);

  // cout << "\nORIGINAL\n\n\n\n";

  // //test_graph.display();

  // cost = test_graph.get_cost();
  // cout << "\n\nCOST--> summ: " << cost.first << ", corr: "<< cost.second << ", TOTAL: "<< cost.first + cost.second <<" \n\n";

  // vector<int> mrg_ids{1,3,4,5};
  //  //vector<int> mrg_ids{50,51,52,53,54,55,56,57,58,59,60};

  // int merged_id = test_graph.merge_snodes(mrg_ids);

  // vector<int> mrg_ids2{7,8};
  //  //vector<int> mrg_ids{50,51,52,53,54,55,56,57,58,59,60};

  // int merged_id2 = test_graph.merge_snodes(mrg_ids2);

  // cout << "\n\n\n\nMERGED\n\n\n\n";

  // //test_graph.display();

  // cost = test_graph.get_cost();
  // cout << "\n\nCOST--> summ: " << cost.first << ", corr: "<< cost.second << ", TOTAL: "<< cost.first + cost.second<<" \n\n";

  // cout << "\nEXPANDING\n\n\n\n";

  // vector<int> expand = {merged_id, merged_id2};
  // test_graph.expand_snodes( expand);

  // //test_graph.display();

  // cost = test_graph.get_cost();
  // cout << "\n\nCOST--> summ: " << cost.first << ", corr: "<< cost.second << ", TOTAL: "<< cost.first + cost.second<<" \n\n";

 /////////////////////////////////////////////////////////////////////////////////////////////
                    /*  TEST HASHING   */
 ////////////////////////////////////////////////////////////////////////////////////////////


  // graph test_graph = graph(opts.graph, opts.types, opts.no_mult);

  // auto fwd_graph = test_graph.get_fwd_graph();
  // auto rev_graph = test_graph.get_rev_graph();

  // cout << "\n\nHashind forward graph..\n\n";
  // auto fwd_list = increment_lsh(fwd_graph, opts.rows, opts.bands );
  // cout << "\n\nHashind reverse graph..\n\n";  
  // auto rev_list = increment_lsh(rev_graph, opts.rows, opts.bands );
  // cout << "\nMerging lists and generating candidates\n\n";
  // auto cand = get_candidates(fwd_list, rev_list, fwd_graph, rev_graph, opts.K);

  // fwd_graph.clear();
  // rev_graph.clear();

/////////////////////////////////////////////////////////////////////////////////////////////
                    /* TEST MERGING  */
/////////////////////////////////////////////////////////////////////////////////////////////

  // pair<double,double> cost;

  // graph test_graph = graph(opts.graph, opts.types, opts.no_mult);

  // cout << "\nORIGINAL\n\n\n\n";

  // //test_graph.display();

  // cost = test_graph.get_cost();
  // cout << "\n\nCOST--> summ: " << cost.first << ", corr: "<< cost.second << ", TOTAL: "<< cost.first + cost.second <<" \n\n";

  // // vector<int> mrg_ids{50,51,52,53,54,55,56,57,58,59,60};
  // vector<int> mrg_ids{2,3,4,5};

  // test_graph.merge_snodes(mrg_ids);

  // cout << "\n\n\n\nMERGED\n\n\n\n";

  // //test_graph.display();

  // cost = test_graph.get_cost();
  // cout << "\n\nCOST--> summ: " << cost.first << ", corr: "<< cost.second << ", TOTAL: "<< cost.first + cost.second<<" \n\n";


/////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////

  // pair<double,double> cost;

  // graph test_graph = graph(opts.graph, opts.types, opts.no_mult);

  // cout << "\nORIGINAL\n\n\n\n";

  // test_graph.display();

  // cost = test_graph.get_cost();
  // cout << "\n\nCOST--> summ: " << cost.first << ", corr: "<< cost.second << ", TOTAL: "<< cost.first + cost.second <<" \n\n";

  // vector<int> mrg_ids{1,2,3,5};

  // test_graph.merge_snodes(mrg_ids);

  // cout << "\n\n\n\nMERGED\n\n\n\n";

  // test_graph.display();

  // cost = test_graph.get_cost();
  // cout << "\n\nCOST--> summ: " << cost.first << ", corr: "<< cost.second << ", TOTAL: "<< cost.first + cost.second<<" \n\n";

  // vector<int> mrg_ids2{10,4};

  // test_graph.merge_snodes(mrg_ids2);

  // cout << "\n\n\n\nMERGED-2\n\n\n\n";

  // test_graph.display();

  // cost = test_graph.get_cost();
  // cout << "\n\nCOST--> summ: " << cost.first << ", corr: "<< cost.second << ", TOTAL: "<< cost.first + cost.second <<" \n\n";

  // vector<int> mrg_ids3{6,7,8,9};

  // test_graph.merge_snodes(mrg_ids3);

  // cout << "\n\n\n\nMERGED-3\n\n\n\n";

  // test_graph.display();

  // cost = test_graph.get_cost();
  // cout << "\n\nCOST--> summ: " << cost.first << ", corr: "<< cost.second << ", TOTAL: "<< cost.first + cost.second <<" \n\n";

  // vector<int> mrg_ids4{11, 12};

  // test_graph.merge_snodes(mrg_ids4);

  // cout << "\n\n\n\nMERGED-4\n\n\n\n";

  // test_graph.display();

  // cost = test_graph.get_cost();
  // cout << "\n\nCOST--> summ: " << cost.first << ", corr: "<< cost.second << ", TOTAL: "<< cost.first + cost.second <<" \n\n";

  return 0;
}
