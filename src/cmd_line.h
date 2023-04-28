#ifndef __CMD_LINE_H
#define __CMD_LINE_H

/// command line information and global parameters
class cmd_line {
public:
    
    /// terms to bound the optimization
    double t_max ;
    double t_min ;
    
    /// to bound proportion search
    double p_max ;
    double p_min ;
    
    /// to create intial simplex points
    double t_length ;
    double p_length ;
    
    /// number of restarts
    int n_restarts ;
    
    /// proportion ancestry for 0-n must sum to 1
    /// these are therefore the final ancestry proportion, not necessary the proportion that fluxed if additional pulses occured closer to the present
    vector<double> ancestry_proportion ;
    
    /// store relevant ancestry information
    vector<pulse> ancestry_pulses ;
    
    /// diploid effective population size ( i.e. 2n )
    double ne ;
    
    /// tolerance for parameter search
    double tolerance ;
    
    /// minimum recombinational distance between markers
    double minimum_distance ;
    
    /// error rates for reads (if read based) or genotypes (if genotype based)
    double error_rate ;
    
    /// bool sample is expressed as genotypes, not read counts
    bool genotype ;
    
    /// ancestral genotype frequencies are fixed
    bool ancestral_fixed ;
    
    /// viterbi output
    /// caution: not recommended for more samples of ploidy > 1
    bool viterbi ;
    
    /// number of digits of precision to include
    int precision ;
    
    /// output actual pulses rather than ancestry states
    bool output_pulses ;
    
    /// error rates specifed
    bool error_rates ; 
    
    /// input file name
    string input_file ;
    
    /// sample file
    string sample_file ;
    
    /// bootstrap
    int n_bootstraps ;
    int block_size ; 


    /// model file name, defines MLS models
    string model_file;

    /// Multi-Locus selections or model fits
    vector<Search> mls_searches;

    /// Threads used in MLS fitting
    int cores = 1;

    /// Output more info to command line when fitting MLS models
    bool verbose_stderr = false;

    /// Radius around which sites are accounted for in MLS fast lnl calculation
    double fast_transitions_radius_in_morgans = 0.02;

    /// range of lnL in simplex after which first stage of searching an MLS model is stopped
    double search_stage_1_threshold = 5;

    /// range of lnL in simplex after which second stage of searching an MLS model is stopped
    double search_stage_2_threshold = 1;

    /// Adjacent pairs of sampled sites skipped in fast lnl calculation of MLS model
    int sampled_pair_skips = 4;

    // Output relative fitnesses of genotype rather than selection coefficients
    bool output_relative_fitnesses = false;
    
    /// read relevant information
    void read_cmd_line ( int argc, char *argv[] ) ;

} ;

#endif

