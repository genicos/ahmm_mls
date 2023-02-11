#ifndef SELECTION_OPTIMIZATION_CONTEXT
#define SELECTION_OPTIMIZATION_CONTEXT
#include <vector>

class selection_opt{
public:
    vector<double> n_recombs;
    cmd_line options;
    //nelder_mead optimizer;
    vector<markov_chain> markov_chain_information ;
    map<int, vector<vector< map< vector<transition_information>, double > > > > transition_matrix_information ;
    vector<int> position ;
    vector<double> morgan_position;
    double chrom_size;

    double neutral_lnl;
    double fast_neutral_lnl;

    map<int,vector<vector<int> > > state_list ;


    vector<double> restricted_search_sites;
    vector<string> chromosomes ;


    selection_opt(vector<double> neutral_recombs, cmd_line o, vector<markov_chain> mci, map<int, vector<vector<map<vector<transition_information>,double>>>> tmi, vector<int> pos, map<int,vector<vector<int>>> &states, vector<string> chroms){
        n_recombs = neutral_recombs;
        options = o;
        markov_chain_information = mci;
        transition_matrix_information = tmi;
        position = pos;
        state_list = state_list;
        chromosomes = chroms;
        
        morgan_position.resize(n_recombs.size());
        double sum = 0;
        for(uint i = 0; i < n_recombs.size(); i++){
            sum += n_recombs[i];
            morgan_position[i] = sum;
        }

        chrom_size = sum;

    }
    selection_opt(){}

    void examine_models();

    void set_context();
    
};

#endif