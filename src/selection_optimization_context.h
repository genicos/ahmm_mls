#ifndef SELECTION_OPTIMIZATION_CONTEXT
#define SELECTION_OPTIMIZATION_CONTEXT
#include <vector>


// Object containing info neccessary to perform MLS searches and calculate lnL
class selection_opt {
public:

    // Recombination rate between each adjacent pair of sampled sites
    vector<double> n_recombs;

    // Command line arguments
    cmd_line options;
    
    // Markov chain info for lnl calculation
    vector<markov_chain> markov_chain_information ;
    map<int, vector<vector< map< vector<transition_information>, double > > > > transition_matrix_information ;

    // Position of each sampled site in bp
    vector<int> position ;

    // Position of each sampled site in morgans
    vector<double> morgan_position;

    double chrom_size;

    // State list for combining probabilities
    map<int,vector<vector<int> > > state_list ;

    // Name of sampled chromosomes
    vector<string> chromosomes ;

    // Current MLS search
    Search curr_search;

    selection_opt(vector<double> &neutral_recombs, cmd_line &o, vector<markov_chain> &mci, map<int, vector<vector<map<vector<transition_information>,double>>>> &tmi, vector<int> &pos, map<int,vector<vector<int>>> &states, vector<string> &chroms){
        n_recombs = neutral_recombs;
        options = o;
        markov_chain_information = mci;
        transition_matrix_information = tmi;
        position = pos;
        chromosomes = chroms;
        
        // Calculate morgan positions of each sampled site
        morgan_position.resize(n_recombs.size());
        double sum = 0;
        for(uint i = 0; i < n_recombs.size(); i++){
            sum += n_recombs[i];
            morgan_position[i] = sum;
        }

        chrom_size = sum;

    }
    selection_opt(){}

    // Examing MLS models
    void examine_models();
    
};

#endif