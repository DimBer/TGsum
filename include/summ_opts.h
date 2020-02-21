#ifndef SUMM_OPTS_H
#define SUMM_OPTS_H

#include <boost/program_options.hpp>

#define DEFAULT_GRAPH "data/test_graph"
#define DEFAULT_TYPES "data/test_types"
#define DEFAULT_SAVE_TO "IGNORE"
#define DEFAULT_SAVE_ALL_TO "IGNORE"
#define DEFAULT_NTHREADS 4
#define DEFAULT_BANDS 10
#define DEFAULT_NITERS 1
#define DEFAULT_ROWS 5 
#define DEFAULT_K 20

struct my_opt{
    std::string graph;
    std::string types;
    size_t bands;
    size_t rows;
    size_t K;
    size_t niters;
    int nthreads;
    bool no_mult;
    bool no_type;
    bool more_greedy;
    bool more_greedy_light;
    bool save_cands;
    bool wis_based;
    std::string save_to;
    std::string save_all_to;
};

namespace po = boost::program_options;

struct my_opt get_options(int argc, char ** argv)
{
    
    struct my_opt opts;
    opts.no_mult = false;
    opts.no_type = false;
    opts.more_greedy = false;
    opts.save_cands = false;
    opts.more_greedy_light = false;
    opts.wis_based = false;

    //Define options
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("graph", po::value<std::string>()->default_value(DEFAULT_GRAPH), "Set path to edge-list input graph file")
        ("types", po::value<std::string>()->default_value(DEFAULT_TYPES), "Set path to input node-type file")
        ("nthreads",po::value<int>()->default_value(DEFAULT_NTHREADS),"Set number of threads")
        ("bands",po::value<size_t>()->default_value(DEFAULT_BANDS),"Set number of lsh bands")
        ("rows",po::value<size_t>()->default_value(DEFAULT_ROWS),"Set number of lsh rows")
        ("K",po::value<size_t>()->default_value(DEFAULT_K),"Maximum number of candidate sets")        
        ("no-mult", po::bool_switch(&opts.no_mult), "Ignore multiplicities.")
        ("no-type", po::bool_switch(&opts.no_type), "Ignore types.")
        ("save-cand", po::bool_switch(&opts.save_cands), "Save candidate sets to file.")        
        ("more-greed", po::bool_switch(&opts.more_greedy), "Use the more greedy search method.")    
        ("more-greed-light",po::bool_switch(&opts.more_greedy_light), "Use the graph-assisted fast greedy search method.")
        ("wis-based",po::bool_switch(&opts.wis_based), "Use the weight-independent-set (WIS) approach.")        
        ("save-to",po::value<std::string>()->default_value(DEFAULT_SAVE_TO), "File path to store summary graph at the end." )
        ("save-all-to",po::value<std::string>()->default_value(DEFAULT_SAVE_ALL_TO), "Directory path to store summary graph after every succesfull merging." )
        ("niters",po::value<size_t>()->default_value(DEFAULT_NITERS),"Number or times summarization is repeated with different seeds." )
    ;

    // Parse command line
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // Convert to my options
    opts.graph = vm["graph"].as<std::string>();
    opts.types = vm["types"].as<std::string>();
    opts.nthreads = vm["nthreads"].as<int>();
    opts.rows = vm["rows"].as<size_t>();
    opts.bands = vm["bands"].as<size_t>();
    opts.K = vm["K"].as<size_t>();
    opts.save_to = vm["save-to"].as<std::string>();
    opts.save_all_to = vm["save-all-to"].as<std::string>();
    opts.niters = vm["niters"].as<size_t>();

    return opts;
}

#endif