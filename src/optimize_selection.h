#ifndef OPTIMIZE_SELECTION
#define OPTIMIZE_SELECTION
#include <vector>
#include <math.h>

enum parameter_type {admix_frac, time_since_admix, location, dominance_coeff, selection_coeff};

#include <time.h>

double get_wall_time() {
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        // error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


selection_opt context;

void selection_opt::set_context() {
    context = *this;
}

//TODO move these functions into selection_opt so i dont have to deal with all this global nonsense
//TODO remove all talk of m and t searching

bool searching_time = false;
bool searching_m    = false;

Search global_search;


vector<vector<vector<double>>> get_local_genotypes (vector<mat> model) {
    
    map<int,vector<mat> > transition_matrix ;

    
    for ( int m = 0 ; m < context.markov_chain_information.size() ; m ++ ) {
        
        alt_create_transition_matrix( transition_matrix, context.transition_matrix_information[context.markov_chain_information.at(m).number_chromosomes], context.n_recombs, context.position, context.markov_chain_information.at(m).number_chromosomes, model ) ;
            
        for ( int p = 0 ; p < context.markov_chain_information[m].ploidy_switch.size() ; p ++ ) {
            alt_create_transition_matrix( transition_matrix, context.transition_matrix_information[context.markov_chain_information[m].ploidy_switch[p]], context.n_recombs,  context.position, context.markov_chain_information[m].ploidy_switch[p], model ) ;
        }
    }
    
    vector<mat> interploidy_transitions;


    int sample_count = context.markov_chain_information.size();
    int site_count = context.markov_chain_information[0].alphas.size();

    
    //populating alphas
    for ( int m = 0 ; m < context.markov_chain_information.size() ; m ++ ) {
        context.markov_chain_information[m].compute_forward_probabilities(  transition_matrix, interploidy_transitions) ;
    }
    
    //populating betas
    for ( int m = 0 ; m < context.markov_chain_information.size() ; m ++ ) {
        context.markov_chain_information[m].compute_backward_probabilities( transition_matrix, interploidy_transitions ) ;
    }

    //vector<double> data_ancestry(context.markov_chain_information[0].alphas.size());
    vector<vector<vector<double>>> genotype_posteriors(context.markov_chain_information.size());

    
    //looping through samples
    for ( int m = 0 ; m < context.markov_chain_information.size() ; m ++ ) {

        vector<vector<double>> sample_ancestry(context.markov_chain_information[0].alphas.size());

        //looping through sites
        for( int i = 0; i < sample_ancestry.size(); i++){

            vec smoothed_probs = context.markov_chain_information[m].alphas[i] % context.markov_chain_information[m].betas[i] ;
            normalize( smoothed_probs ) ;

            //int ploidy = statecount2ploidy[smoothed_probs.size()];

            //ploidy 1
            if(smoothed_probs.size() == 2) {
                sample_ancestry[i] = conv_to< vector<double> >::from(smoothed_probs);
                //sample_ancestry[i] += smoothed_probs[0];
            }

            //ploidy 2
            if(smoothed_probs.size() == 3) {
                sample_ancestry[i] = conv_to< vector<double> >::from(smoothed_probs);
                //sample_ancestry[i] += smoothed_probs[0];
                //sample_ancestry[i] += smoothed_probs[1] * 0.5;
            }

            //TODO generalize ploidy here
            /*
            data_ancestry[i] += smoothed_probs[0];
            for(int i = 1; i < smoothed_probs.size() - 1; i++) {
                data_ancestry[i] += smoothed_probs[i] * (smoothed_probs.size() - i)/(smoothed_probs.size() - 1);
            }
            */
        }

        genotype_posteriors[m] = sample_ancestry;

    }
    
    return genotype_posteriors;
}








vector<double> get_local_ancestry (vector<mat> model) {
    
    vector<vector<vector<double>>> genotypes = get_local_genotypes(model);

    vector<double> data_ancestry(context.markov_chain_information[0].alphas.size());
    
    //looping through samples
    for ( int m = 0 ; m < context.markov_chain_information.size() ; m ++ ) {

        //looping through sites
        for( int i = 0; i < data_ancestry.size(); i++){

            //ploidy 1
            if(genotypes[m][i].size() == 2) {
                data_ancestry[i] += genotypes[m][i][0];
            }

            //ploidy 2
            if(genotypes[m][i].size() == 3) {
                data_ancestry[i] += genotypes[m][i][0];
                data_ancestry[i] += genotypes[m][i][1] * 0.5;
            }

        }

    }
    
    for( int i = 0; i < data_ancestry.size(); i++){
        data_ancestry[i] /= context.markov_chain_information.size();
    }
    
    return data_ancestry;
}




void prepare_selection_info(vector<double> &parameters, vector<double> &selection_recomb_rates, vector<vector<double>> &fitnesses){

    int selected_sites_count = parameters.size() / 3;
    selection_recomb_rates.resize(selected_sites_count + 1);
    fitnesses.resize(selected_sites_count);
    
    if(parameters.size() > 0){
        

        if(context.options.verbose_stderr) {
            cerr << "\nTesting parameters:\n";
            if (global_search.search_m) {
                cerr << "m: " << parameters[0] << "\n";
            }
            if (global_search.search_t) {
                cerr << "time:    " << parameters[global_search.search_m] << "\n";
            }
            for (uint i = 0; i < selected_sites_count; i++) {
                cerr << "selection site: " << parameters[3*i + 0 + global_search.search_m + global_search.search_t] << " with fitness: " << parameters[3*i + 1 + global_search.search_m + global_search.search_t] << ",1," << parameters[3*i + 2 + global_search.search_m + global_search.search_t] << "\n";
            }
        }else{
            if (searching_time) {
                cerr << "m: " << parameters[0] << "\n";
            }
            if (searching_m) {
                cerr << "time:    " << parameters[global_search.search_m] << "\n";
            }
            for (uint i = 0; i < selected_sites_count; i++) {
                cerr << parameters[3*i + 0 + global_search.search_m + global_search.search_t] << "\t" << parameters[3*i + 1 + global_search.search_m + global_search.search_t] << "\t" << parameters[3*i + 2 + global_search.search_m + global_search.search_t] << "\n";
            }
        }

        //Sort selected sites and fitnesses///////////////////////////
        struct Selected_pair{                                       //
            double site;
            vector<double> fitness;

            bool operator<(const Selected_pair x) const
                { return site < x.site;}
        };
        
        //Fill vector of structs
        vector<Selected_pair> selected_pairs(selected_sites_count);
        for(int i = 0; i < selected_sites_count; i++){
            Selected_pair ss;
            ss.site = parameters[i*3 + global_search.search_m + global_search.search_t];
            ss.fitness.resize(3);
            ss.fitness[0] = parameters[i*3 + 1 + global_search.search_m + global_search.search_t];
            ss.fitness[1] = 1;
            ss.fitness[2] = parameters[i*3 + 2 + global_search.search_m + global_search.search_t];

            selected_pairs[i] = ss;
        }
              
        //sort
        sort(selected_pairs.begin(), selected_pairs.end());         //
        //////////////////////////////////////////////////////////////

        selection_recomb_rates[0] = 1;

        double last = 0;
        for (uint i = 0; i < selected_pairs.size(); i++){
            selection_recomb_rates[i] = selected_pairs[i].site - last;
            last = selected_pairs[i].site;
        }

        selection_recomb_rates[selected_pairs.size()] = 1 - last;

        for(uint i = 0; i < selected_sites_count; i++){
            fitnesses[i] = selected_pairs[i].fitness;
        }
        
    }else{
        selection_recomb_rates[0] = 1;
    }

}




double compute_lnl(vector<mat> &transition_matrices){

    map<int,vector<mat> > transition_matrix ;
    
    // Compute transition matricies for different ploidies
    for ( int m = 0 ; m < context.markov_chain_information.size() ; m ++ ) {
        
        alt_create_transition_matrix( transition_matrix, context.transition_matrix_information[context.markov_chain_information.at(m).number_chromosomes], context.n_recombs, context.position, context.markov_chain_information.at(m).number_chromosomes, transition_matrices ) ;
        
        for ( int p = 0 ; p < context.markov_chain_information[m].ploidy_switch.size() ; p ++ ) {
            alt_create_transition_matrix( transition_matrix, context.transition_matrix_information[context.markov_chain_information[m].ploidy_switch[p]], context.n_recombs, context.position, context.markov_chain_information[m].ploidy_switch[p], transition_matrices ) ;
        }

    }
    
    
    vector<mat> interploidy_transitions;

    
    double lnl = 0 ;
    
    // Sum up log likelihoods for each panel
    for ( int m = 0 ; m < context.markov_chain_information.size() ; m ++ ) {
        
        lnl += context.markov_chain_information[m].compute_lnl( transition_matrix, interploidy_transitions) ;
    }
    
    return lnl;
}





vector<mat> last_calculated_transition_matricies;

double to_be_optimized (vector<double> parameters) {

    double timer = get_wall_time();
    
    vector<double> selection_recomb_rates;
    vector<vector<double>> fitnesses;
    

    prepare_selection_info(parameters, selection_recomb_rates, fitnesses);

    int cores = context.options.cores;

    double m        = (global_search.search_m) ?      parameters[0] : global_search.start_m;
    int generations = (global_search.search_t) ? (int)parameters[global_search.search_m] : global_search.start_t;

    vector<mat> transition_matrices = calculate_transition_rates (
        context.n_recombs,
        selection_recomb_rates,
        fitnesses,
        m,
        generations,
        cores
    );
    
    last_calculated_transition_matricies = transition_matrices;
    
    double lnl = compute_lnl(transition_matrices);

    if(context.options.verbose_stderr){
        cerr << "lnl = " << setprecision(15) << lnl << "\n";
        cerr << "TIME PASSED IN ONE ITERATION: " << (get_wall_time() - timer) << "\n";
    }else{
        cerr << "lnl\t" << setprecision(15) << lnl << "\n";
    }
    
    return lnl;
}





vector<mat> neutral_transition_matrices;

double to_be_optimized_only_near_sites(vector<double> parameters) {
    
    double timer = get_wall_time();
    
    vector<double> selection_recomb_rates;
    vector<vector<double>> fitnesses;

    prepare_selection_info(parameters, selection_recomb_rates, fitnesses);
    
    int cores = context.options.cores;
    

    if(context.options.use_model_file) // TODO, chnage how cores affects model files
        cores = 1;

    
    double m        = (global_search.search_m) ?      parameters[0] : global_search.start_m;
    int generations = (global_search.search_t) ? (int)parameters[global_search.search_m] : global_search.start_t;
    
    vector<mat> transition_matrices = alternative_fast_transition_rates (
        context.n_recombs,
        selection_recomb_rates,
        fitnesses,
        m,
        generations,
        cores,
        context.options.fast_transitions_radius_in_morgans
    );

    double lnl = compute_lnl(transition_matrices);
    

    if(context.options.verbose_stderr) {
        cerr << "lnl ratio = " << setprecision(15) << lnl << "\n";
        cerr << "TIME PASSED IN ONE ITERATION: " << (get_wall_time() - timer) << "\n";
    }else{
        cerr << "lnl ratio\t" << setprecision(15) << lnl << "\n";
    }
    
    return lnl;
}













void setup_searches(Search this_search) {
    global_search = this_search;
}

vector<double> convert_parameters_to_long_form(vector<double> parameters) {
    vector<double> new_params(0);

    int i = 0;
    if (global_search.search_m) {
        new_params.push_back(parameters[i++]);
    }

    if(global_search.search_t) {
        new_params.push_back(parameters[i++]);
    }

    for(int j = 0; j < global_search.search_l.size(); j++) {

        if(global_search.search_l[j]) {
            new_params.push_back(parameters[i++]);
        }else{
            new_params.push_back(global_search.start_l[j]);
        }

        bool h_and_s = true; // Not all fitnesses can be expressed in terms of h and s
        double h = 0;
        double s = 0;

        if(global_search.search_h[j]){
            h = parameters[i++];
        }else{
            h = global_search.start_h[j];
        }

        if(global_search.search_s[j]){
            s = parameters[i++];
        }else{
            s = global_search.start_s[j];
        }

        if(h_and_s){
            new_params.push_back(1 / (1 - h * s));
            new_params.push_back((1 - s) / (1 - h * s));
        }
    }


    return new_params;
}

double general_to_be_optimized(vector<double> parameters){
    parameters = convert_parameters_to_long_form(parameters);

    return to_be_optimized(parameters);
}

double general_to_be_optimized_fast(vector<double> parameters){
    parameters = convert_parameters_to_long_form(parameters);

    return to_be_optimized_only_near_sites(parameters);
}











vector<double> multi_level_optimization(
    double chrom_size,
    nelder_mead &opt,
    vector<double> &given_parameters,
    vector<vector<vector<double>>> &bottle_necks,
    double (*to_be_optimized_function) (vector<double>),
    vector<parameter_type> parameter_types
){
    vector<double> best_parameters;
    double best_ratio = -DBL_MAX;

    vector<double> found_parameters;
    

    for(int j = 0; j < bottle_necks.size(); j++){
        for(int k = 0; k < bottle_necks[j].size(); k++){
            for(int l = 0; l < bottle_necks[j][k][0]; l++){
                
                cerr << "\n SEARCH " << j + 1 << "/" << bottle_necks.size() << " " << k + 1 << "/" << bottle_necks[j].size() << " " << l + 1 << "/" << bottle_necks[j][k][0] << "\n";
                
                

                vector<double> center_point = given_parameters;
                vector<double> scales(given_parameters.size());


                for(int i = 0; i < given_parameters.size(); i++) {
                    switch(parameter_types[i]){
                        case admix_frac:
                            scales[i] = bottle_necks[j][k][2];
                            break;
                        case time_since_admix:
                            scales[i] = bottle_necks[j][k][1] * 1000;
                            break;
                        case location:
                            scales[i] = bottle_necks[j][k][1];
                            break;
                        case dominance_coeff:
                            scales[i] = bottle_necks[j][k][2] * 5;
                            break;
                        case selection_coeff:
                            scales[i] = bottle_necks[j][k][2];
                            break;
                    }
                }


                opt.populate_points(given_parameters.size(), 1, center_point, scales);

                

                //random reflections
                for(uint i = 0; i < center_point.size(); i++){
                    if(rand() < RAND_MAX/2){
                        for(uint j = 0; j < opt.points.size(); j++){
                            opt.points[j][i] = 2*center_point[i] - opt.points[j][i];
                        }
                    }
                }
                

                opt.enforce_bounds();



                opt.calculate_points(to_be_optimized_function);

                double initial_range = opt.max_value - opt.min_value;

                double stopping_point = min(initial_range / 20 * bottle_necks[j][k][3], bottle_necks[j][k][3]);

                while((opt.max_value - opt.min_value) > stopping_point && opt.repeated_shrinkages < 4){
                    opt.iterate(to_be_optimized_function);

                    if(context.options.verbose_stderr){
                        cerr << "nelder-mead simplex_size: " << opt.simplex_size() << " output range: " << opt.max_value - opt.min_value << "\n";
                    }
                }

                found_parameters = opt.points[opt.max_index];

                if(opt.max_value > best_ratio) {
                    best_ratio = opt.max_value;
                    best_parameters = found_parameters;
                }
            }
        }

        given_parameters = best_parameters;

        best_ratio = -DBL_MAX;
    }

    return best_parameters;
}


#endif