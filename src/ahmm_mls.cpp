/*

 copyright: Russ Corbett-Detig
            rucorbet@ucsc.edu

 MLS extension created by:  Nicolas Ayala
                            nmayala@ucsc.edu
 
 This is software distributed under the gnu public license version 3.
 
*/

/// headers
#include <iostream>
#include <vector>
#include <map>
#include <time.h>
#include <string>
#include <fstream>
#include <algorithm>
using namespace std ;

/// linear algebra library is armadillo
#define ARMA_NO_DEBUG
#include <armadillo>
using namespace arma ;

/// our header files in /src directory
#include "print_usage.h"
#include "factorial.h"
#include "nchoosek.h" 
#include "subsample.h" 
#include "multichoose.h"
#include "multipermute.h"
#include "normalize.h"
#include "ancestry_pulse.h"
#include "ploidy_path.h" 
#include "markov_chain.h"
#include "read_samples.h" 
#include "pulses_to_ancestry.h" 
#include "compute_forward.h"
#include "compute_backward.h"
#include "forward_backward.h"
#include "viterbi.h" 
#include "transition_information.h"
#include "exponentiate_matrix.h"
#include "cmd_line.h"
#include "create_transition_rates.h"
#include "read_cmd_line.h"
#include "evaluate_vertex.h"
#include "check_vertex.h"
#include "sort_vertices.h"
#include "create_pulses.h" 
#include "create_states.h"
#include "input_line.h"
#include "distribute_alleles.h" 
#include "binomial.h"
#include "read_emissions.h"
#include "genotype_emissions.h"
#include "read_input.h"
#include "nelder_mead.h"
#include "golden_search.h"
#include "bootstrap.h" 
#include <iomanip>

#include "selection_optimization_context.h"
#include "simulate_selective_admixture.h"
#include "new_read_input.h"
#include "nelder_mead_general.h"
#include "optimize_selection.h"



#include "read_sites.h"
#include "grid_search.h"
#include "site_examining.h"

#include "read_model_file.h"
#include "model_examining.h"

#include "uninformed_inference.h"


int main ( int argc, char *argv[] ) {
    
    /// time tracking
    clock_t t = clock() ;
    clock_t total = clock() ;
    
    /// seed prng
    srand (t) ;
    
	// read cmd line 
	cmd_line options ;
    cerr << "reading command line" ; t = clock();
	options.read_cmd_line( argc, argv ) ;

    /// chain objects for each sample
    vector<markov_chain> markov_chain_information ;
    
    /// get sample ids and ploidy from input file
    cerr << "\t\t\t\t" << (double) (clock() - t) << " ms\n" << "reading sample ids and ploidy" ; t = clock();
    read_samples( markov_chain_information, options.sample_file, options.viterbi ) ;

    /// create states matrix
    cerr << "\t\t\t" << (double) (clock() - t) << " ms\n" << "creating states matrix" ; t = clock();
    /// store all possible state space arranged by ploidy and then vector of state counts
    map<int,vector<vector<int> > > state_list ;
    /// now create initial state list
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        for ( int p = 0 ; p < markov_chain_information[m].sample_ploidy_path.size() ; p ++ ) {
            create_initial_states( markov_chain_information.at(m).sample_ploidy_path[p].ploidy, options.ancestry_pulses, state_list ) ;
        }
    }
    
	/// read in panels and update matrices
    cerr << "\t\t\t\t" << (double) (clock() - t) << " ms\n" << "reading data and creating emissions matrices\t" ; t = clock() ;
    /// store recombination rates and positions
    vector<int> position ;
    vector<double> recombination_rate ;
    vector<string> chromosomes ;


    
    
    selection_read_file( options, markov_chain_information, state_list, position, recombination_rate, chromosomes ) ;
    
    
    
    recombination_rate[0] = 0;

    
    
    
        
    /// create basic transition information
    cerr << (double) (clock() - t) << " ms" << endl << "computing transition routes\t\t\t" ; t = clock() ;
    /// 3d map to look up by ploidy, start state, end state, and then relevant transition information
    map<int, vector<vector< map< vector<transition_information>, double > > > > transition_matrix_information ;
    for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
        for ( int p = 0 ; p < markov_chain_information[m].sample_ploidy_path.size() ; p ++ ) {
            create_transition_information( markov_chain_information.at(m).sample_ploidy_path[p].ploidy, transition_matrix_information, state_list[markov_chain_information.at(m).sample_ploidy_path[p].ploidy] ) ;
        }
    }
    
    /// create admixture model(s)
    cerr << (double) (clock() - t) << " ms" << endl << "creating initial admixture model(s)\t\t" ; t = clock();
    vector<vector<pulse> > vertices ;
    int nparams = create_pulses( vertices, options ) ;
    
    /// set number of restarts if unspecified default is factorial * 2
    
    if ( options.n_restarts < 0 ) {
        options.n_restarts = factorial[nparams] * 2 ;
    }
    
    
    if(options.uninformed_inference){
        selection_opt selection_optimizer(recombination_rate, options, markov_chain_information, transition_matrix_information, position);

        selection_optimizer.uninformed_inference();
    }

    if(options.use_model_file){
        read_model_file(options);
        selection_opt selection_optimizer(recombination_rate, options, markov_chain_information, transition_matrix_information, position);
            
        selection_optimizer.test_models();
    }
    
    
    if(options.use_site_file){
        read_site_file(options, recombination_rate, position);

        selection_opt selection_optimizer(recombination_rate, options, markov_chain_information, transition_matrix_information, position);
        
        selection_optimizer.examine_sites();
    }


    

    
    
        
    cerr << (double) (clock() - t) << " ms" << endl ;
    cerr << "total run time:\t\t\t" << (double) (clock() - total) << " ms" << endl ;

	return 0 ; 
}
	

