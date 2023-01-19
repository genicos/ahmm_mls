#ifndef OPTIMIZE_SELECTION
#define OPTIMIZE_SELECTION
#include <vector>
#include <math.h>

enum parameter_type {location, selection_coeff, dominance_coeff, time_since_admix};

#include <time.h>
#include <sys/time.h>
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

bool searching_time = false;
bool searching_m    = false;


vector<double> get_local_ancestry (vector<mat> neutral_model) {
    
    map<int,vector<mat> > transition_matrix ;

    
    for ( int m = 0 ; m < context.markov_chain_information.size() ; m ++ ) {
        
        alt_create_transition_matrix( transition_matrix, context.transition_matrix_information[context.markov_chain_information.at(m).number_chromosomes], context.n_recombs, context.position, context.markov_chain_information.at(m).number_chromosomes, neutral_model ) ;
            
        for ( int p = 0 ; p < context.markov_chain_information[m].ploidy_switch.size() ; p ++ ) {
            alt_create_transition_matrix( transition_matrix, context.transition_matrix_information[context.markov_chain_information[m].ploidy_switch[p]], context.n_recombs,  context.position, context.markov_chain_information[m].ploidy_switch[p], neutral_model ) ;
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

    vector<double> data_ancestry(context.markov_chain_information[0].alphas.size());

    
    //looping through samples
    for ( int m = 0 ; m < context.markov_chain_information.size() ; m ++ ) {

        //looping through sites
        for( int i = 0; i < data_ancestry.size(); i++){

            

            vec smoothed_probs = context.markov_chain_information[m].alphas[i] % context.markov_chain_information[m].betas[i] ;
            normalize( smoothed_probs ) ;

            //int ploidy = statecount2ploidy[smoothed_probs.size()];

            //ploidy 1
            if(smoothed_probs.size() == 2) {
                data_ancestry[i] += smoothed_probs[0];
            }

            //ploidy 2
            if(smoothed_probs.size() == 3) {
                data_ancestry[i] += smoothed_probs[0];
                data_ancestry[i] += smoothed_probs[1] * 0.5;
            }


            //TODO generalize ploidy here
            /*
            data_ancestry[i] += smoothed_probs[0];
            for(int i = 1; i < smoothed_probs.size() - 1; i++){
                data_ancestry[i] += smoothed_probs[i] * (smoothed_probs.size() - i)/(smoothed_probs.size() - 1);
            }
            */
        }

    }
    
    for( int i = 0; i < data_ancestry.size(); i++){
        data_ancestry[i] /= sample_count;
    }
    
    return data_ancestry;
}




void prepare_selection_info(vector<double> &parameters, vector<double> &selection_recomb_rates, vector<vector<double>> &fitnesses){

    int selected_sites_count = parameters.size()/3;
    selection_recomb_rates.resize(selected_sites_count + 1);
    fitnesses.resize(selected_sites_count);

    if(selected_sites_count > 0){
        

        if(context.options.verbose_stderr) {
            cerr << "\nTesting parameters:\n";
            if (searching_time){
                cerr << "Time: " << parameters[0] << "\n";
            }
            if (searching_m){
                cerr << "m:    " << parameters[searching_time] << "\n";
            }
            for (uint i = 0; i < selected_sites_count; i++) {
                cerr << "selection site: " << parameters[3*i + 0 + searching_time + searching_m] << " with fitness: " << parameters[3*i + 1 + searching_time + searching_m] << ",1," << parameters[3*i + 2 + searching_time + searching_m] << "\n";
            }
        }else{
            if (searching_time){
                cerr << "Time: " << parameters[0] << "\n";
            }
            if (searching_m){
                cerr << "m:    " << parameters[searching_time] << "\n";
            }
            for (uint i = 0; i < selected_sites_count; i++) {
                cerr << parameters[3*i + 0 + searching_time + searching_m] << "\t" << parameters[3*i + 1 + searching_time + searching_m] << "\t" << parameters[3*i + 2 + searching_time + searching_m] << "\n";
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
            ss.site = parameters[i*3 + searching_time + searching_m];
            ss.fitness.resize(3);
            ss.fitness[0] = parameters[i*3 + 1 + searching_time + searching_m];
            ss.fitness[1] = 1;
            ss.fitness[2] = parameters[i*3 + 2 + searching_time + searching_m];

            selected_pairs[i] = ss;
        }
              
        //sort
        sort(selected_pairs.begin(), selected_pairs.end());         //
        //////////////////////////////////////////////////////////////


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

    if(context.options.use_model_file) // TODO mayhaps delete this
        cores = 1;

    int generations = (searching_time) ? (int)parameters[0] : context.options.generations;
    double m = (searching_m) ? parameters[searching_time] : context.options.m;

    vector<mat> transition_matrices = calculate_transition_rates(
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

    int generations = (searching_time) ? (int)parameters[0] : context.options.generations;
    double m = (searching_m) ? parameters[searching_time] : context.options.m;
    
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
        cerr << "lnl ratio = " << setprecision(15) << lnl - context.fast_neutral_lnl<< "\n";
        cerr << "TIME PASSED IN ONE ITERATION: " << (get_wall_time() - timer) << "\n";
    }else{
        cerr << "lnl ratio\t" << setprecision(15) << lnl - context.fast_neutral_lnl<< "\n";
    }
    
    return lnl;
}





/*
double to_be_optimized_general();
*/







//each site is two parameters,p[0] p[1], which translates to (p[0], 1, p[1])
double to_be_optimized_pop0_dominant(vector<double> parameters){

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]); 
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }
    int start = searching_time + searching_m;

    for(; i < parameters.size(); i++){
        new_parameters.push_back(parameters[i]);
        if((i-start) % 2 == 0){
            new_parameters.push_back(1);
        }
    }

    return to_be_optimized(new_parameters);
}

//each site is two parameters,p[0] p[1], which translates to (p[0], p[1], 1)
double to_be_optimized_pop1_dominant(vector<double> parameters){

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]); 
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }
    int start = searching_time + searching_m;

    for(; i < parameters.size(); i++){
        new_parameters.push_back(parameters[i]);
        if((i-start) % 2 == 1){
            new_parameters.push_back(1);
        }
    }

    return to_be_optimized(new_parameters);
}

//each site is two parameters,p[0] p[1], which translates to (p[0], (1-p[1])/(1-p[1]/2), 1/(1-p[1]/2))
double to_be_optimized_additive(vector<double> parameters) {

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]); 
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }
    int start = searching_time + searching_m;
    
    for(; i < parameters.size(); i+=2) {
        new_parameters.push_back(parameters[i]);
        new_parameters.push_back((1-parameters[i+1])/(1-parameters[i+1]/2));
        new_parameters.push_back(1/(1-parameters[i+1]/2));
    }

    return to_be_optimized(new_parameters);
}


//each site is three parameters,p[0] p[1] p[2], which translates to (p[0], (1-p[1])/(1-p[1]*p[2]), 1/(1-p[1]*p[2]))
double to_be_optimized_h(vector<double> parameters) {

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]); 
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }
    int start = searching_time + searching_m;

    
    for(; i < parameters.size(); i+=3){
        new_parameters.push_back(parameters[i]);
        new_parameters.push_back((1-parameters[i+1])/(1-parameters[i+1]*parameters[i+2]));
        new_parameters.push_back(1/(1-parameters[i+1]*parameters[i+2]));   
    }

    return to_be_optimized(new_parameters);
}







//each site is two parameters,p[0] p[1], which translates to (p[0], 1, p[1])
double to_be_optimized_pop0_dominant_fast(vector<double> parameters){

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]);
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }
    int start = searching_time + searching_m;

    for(; i < parameters.size(); i++){
        new_parameters.push_back(parameters[i]);
        if((i-start) % 2 == 0){
            new_parameters.push_back(1);
        }
    }

    return to_be_optimized_only_near_sites(new_parameters);
}

//each site is two parameters,p[0] p[1], which translates to (p[0], p[1], 1)
double to_be_optimized_pop1_dominant_fast(vector<double> parameters){

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]); 
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }
    int start = searching_time + searching_m;

    for(; i < parameters.size(); i++){
        new_parameters.push_back(parameters[i]);
        if((i-start) % 2 == 1){
            new_parameters.push_back(1);
        }
    }

    return to_be_optimized_only_near_sites(new_parameters);
}


//each site is two parameters,p[0] p[1], which translates to (p[0], (1-p[1])/(1-p[1]/2), 1/(1-p[1]/2))
double to_be_optimized_additive_fast(vector<double> parameters) {

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]); 
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }

    for(; i < parameters.size(); i+=2){
        new_parameters.push_back(parameters[i]);
        new_parameters.push_back((1-parameters[i+1])/(1-parameters[i+1]/2));
        new_parameters.push_back(1/(1-parameters[i+1]/2));
    }

    return to_be_optimized_only_near_sites(new_parameters);
}


//each site is three parameters,p[0] p[1] p[2], which translates to (p[0], (1-p[1])/(1-p[1]*p[2]), 1/(1-p[1]*p[2]))
double to_be_optimized_h_fast(vector<double> parameters) {

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]); 
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }

    for(; i < parameters.size(); i+=3){
        new_parameters.push_back(parameters[i]);
        new_parameters.push_back((1-parameters[i+1])/(1-parameters[i+1]*parameters[i+2]));
        new_parameters.push_back(1/(1-parameters[i+1]*parameters[i+2]));   
    }

    return to_be_optimized_only_near_sites(new_parameters);
}





double to_be_optimized_restricted(vector<double> parameters){

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]); 
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }
    int start = searching_time + searching_m;

    for(; i < parameters.size()/2; i++){
        
        new_parameters.push_back(context.restricted_search_sites[i - start]);
        new_parameters.push_back(parameters[(i-start)*2 + 0 + start]);
        new_parameters.push_back(parameters[(i-start)*2 + 1 + start]);
        
    }

    return to_be_optimized(new_parameters);
}


double to_be_optimized_restricted_fast(vector<double> parameters) {

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]); 
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }
    int start = searching_time + searching_m;

    for(; i < parameters.size()/2; i++){
        
        new_parameters.push_back(context.restricted_search_sites[i - start]);
        new_parameters.push_back(parameters[(i-start)*2 + 0 + start]);
        new_parameters.push_back(parameters[(i-start)*2 + 1 + start]);
        
    }

    return to_be_optimized_only_near_sites(new_parameters);
}












double to_be_optimized_restricted_additive(vector<double> parameters){

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]); 
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }
    int start = searching_time + searching_m;

    for(; i < parameters.size(); i++){
        
        new_parameters.push_back(context.restricted_search_sites[i - start]);
        new_parameters.push_back((1-parameters[i])/(1-parameters[i]/2));
        new_parameters.push_back(1/(1-parameters[i]/2));
        
    }

    return to_be_optimized(new_parameters);
}

//each site is one parameters,p[0] which translates to (NA, (1-p[0])/(1-p[0]/2), 1/(1-p[0]/2))
double to_be_optimized_restricted_additive_fast(vector<double> parameters) {

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]); 
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }
    int start = searching_time + searching_m;

    for(; i < parameters.size(); i++){
        
        new_parameters.push_back(context.restricted_search_sites[i - start]);
        new_parameters.push_back((1-parameters[i])/(1-parameters[i]/2));
        new_parameters.push_back(1/(1-parameters[i]/2));
        
    }

    return to_be_optimized_only_near_sites(new_parameters);
}



double to_be_optimized_restricted_h(vector<double> parameters) {

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]); 
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }
    int start = searching_time + searching_m;
    

    for(; i < parameters.size(); i+=2){
        new_parameters.push_back(context.restricted_search_sites[(i-start)/2]);
        new_parameters.push_back((1-parameters[i])/(1-parameters[i]*parameters[i+1]));
        new_parameters.push_back(1/(1-parameters[i]*parameters[i+1]));   
    }

    return to_be_optimized(new_parameters);
}

double to_be_optimized_restricted_h_fast(vector<double> parameters) {

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]); 
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }
    int start = searching_time + searching_m;

    for(; i < parameters.size(); i+=2){
        new_parameters.push_back(context.restricted_search_sites[(i - start)/2]);
        new_parameters.push_back((1-parameters[i])/(1-parameters[i]*parameters[i+1]));
        new_parameters.push_back(1/(1-parameters[i]*parameters[i+1]));   
    }

    return to_be_optimized_only_near_sites(new_parameters);
}







// (p[0]) => (NA, 1, p[0])
double to_be_optimized_restricted_dom0(vector<double> parameters){

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]); 
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }

    for(; i < parameters.size(); i++){
        
        new_parameters.push_back(context.restricted_search_sites[i]);
        new_parameters.push_back(1);
        new_parameters.push_back(parameters[i]);
        
    }

    return to_be_optimized(new_parameters);
}


double to_be_optimized_restricted_dom0_fast(vector<double> parameters) {

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]); 
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }

    for(; i < parameters.size(); i++){
        
        new_parameters.push_back(context.restricted_search_sites[i]);
        new_parameters.push_back(1);
        new_parameters.push_back(parameters[i]);
        
    }

    return to_be_optimized_only_near_sites(new_parameters);
}





// (p[0]) => (NA, p[0], 1)
double to_be_optimized_restricted_dom1(vector<double> parameters){

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]); 
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }

    for(; i < parameters.size(); i++){
        
        new_parameters.push_back(context.restricted_search_sites[i]);
        new_parameters.push_back(parameters[i]);
        new_parameters.push_back(1);
        
    }

    return to_be_optimized(new_parameters);
}


double to_be_optimized_restricted_dom1_fast(vector<double> parameters) {

    vector<double> new_parameters;

    int i = 0;

    if(searching_time){
        new_parameters.push_back(parameters[0]); 
        i++;
    }
    if(searching_m){
        new_parameters.push_back(parameters[searching_time]); 
        i++;
    }

    for(; i < parameters.size(); i++){
        
        new_parameters.push_back(context.restricted_search_sites[i]);
        new_parameters.push_back(parameters[i]);
        new_parameters.push_back(1);
        
    }

    return to_be_optimized_only_near_sites(new_parameters);
}
































typedef double (*lnl_function)(vector<double> parameters);


lnl_function to_be_optimized_variations (bool fast, bool restricted_site, bool additive, bool dom0, bool dom1, bool search_h, bool time, bool mixf) {

    searching_time = time;
    searching_m    = mixf;

    if(restricted_site){
        if(fast){
            if(additive){
                return &to_be_optimized_restricted_additive_fast;
            }else if(dom0){
                return &to_be_optimized_restricted_dom0_fast;
            }else if(dom1){
                return &to_be_optimized_restricted_dom1_fast;
            }else if(search_h){
                return &to_be_optimized_restricted_h_fast;
            }else{
                return &to_be_optimized_restricted_fast;
            }
        }else{
            if(additive){
                return &to_be_optimized_restricted_additive;
            }else if(dom0){
                return &to_be_optimized_restricted_dom0;
            }else if(dom1){
                return &to_be_optimized_restricted_dom1;
            }else if(search_h){
                return &to_be_optimized_restricted_h;
            }else{
                return &to_be_optimized_restricted;
            }
        }
    }else{
        if(fast){
            if(additive){
                return &to_be_optimized_additive_fast;
            }else if(dom0){
                return &to_be_optimized_pop0_dominant_fast;
            }else if(dom1){
                return &to_be_optimized_pop1_dominant_fast;
            }else if(search_h){
                return &to_be_optimized_h_fast;
            }else{
                return &to_be_optimized_only_near_sites; //TODO i dont like this name
            }
        }else{
            if(additive){
                return &to_be_optimized_additive;
            }else if(dom0){
                return &to_be_optimized_pop0_dominant;
            }else if(dom1){
                return &to_be_optimized_pop1_dominant;
            }else if(search_h){
                return &to_be_optimized_h;
            }else{
                return &to_be_optimized;
            }
        }
    }
}




























//sites to parameters cannot account for time parameter
vector<double> sites_to_parameters(vector<vector<double>> sites, int parameters_per_site = 3){
    if(parameters_per_site == 0){
        parameters_per_site = 3;
    }

    vector<double> parameters(sites.size()*parameters_per_site);

    for(uint i = 0; i < parameters.size(); i++){
        parameters[i] = sites[i/parameters_per_site][i%parameters_per_site];
    }

    return parameters;
}

vector<vector<double>> parameters_to_sites(vector<double> parameters, int parameters_per_site = 3){
    
    vector<vector<double>> sites(parameters.size()/parameters_per_site);
    
    for(int j = 0; j < sites.size(); j++){
        vector<double> site(parameters_per_site);

        for(int i = 0; i < parameters_per_site; i++){
            site[i] = parameters[j*parameters_per_site + i + searching_time + searching_m];
        }
        
        sites[j] = site;
    }
    
    return sites;
}





vector<double> search_site                       (double chrom_size, nelder_mead &opt, vector<double> site, double width, double height, double depth);


vector<double> search_sites                      (double chrom_size, nelder_mead &opt, vector<vector<double>> sites, double width, double height, double depth);
vector<double> search_sites_fast                 (double chrom_size, nelder_mead &opt, vector<vector<double>> sites, double width, double height, double depth);



vector<double> multi_level_optimization(
    double chrom_size,
    nelder_mead &opt,
    vector<double> &given_parameters,
    vector<vector<vector<double>>> &bottle_necks,
    double (*to_be_optimized_function) (vector<double>),
    vector<parameter_type> parameter_types,
    int parameters_per_site = 3
){
    vector<double> best_parameters;
    double best_ratio = -DBL_MAX;

    vector<double> found_parameters;
    

    for(int j = 0; j < bottle_necks.size(); j++){
        for(int k = 0; k < bottle_necks[j].size(); k++){
            for(int l = 0; l < bottle_necks[j][k][0]; l++){
                
                //if(context.options.verbose_stderr){
                    cerr << "\n SEARCH " << j + 1 << "/" << bottle_necks.size() << " " << k + 1 << "/" << bottle_necks[j].size() << " " << l + 1 << "/" << bottle_necks[j][k][0] << "\n";
                //}
                

                vector<double> center_point(given_parameters.size());
                vector<double> scales(given_parameters.size());

                int number_of_sites = (given_parameters.size() - searching_time - searching_m)/parameters_per_site;
                if(searching_time){
                    scales[0] = bottle_necks[j][k][1] * 1000;
                    center_point[0] = given_parameters[0];
                }
                if(searching_m){
                    scales[searching_time] = bottle_necks[j][k][2];
                    center_point[searching_time] = given_parameters[searching_time];
                }

                for(uint i = 0; i < number_of_sites; i++){
                    for(int h = 0; h < parameters_per_site; h++){

                        center_point[i*parameters_per_site + h + searching_time + searching_m] = given_parameters[i*parameters_per_site + h + searching_time + searching_m];

                        if(parameter_types[h] == location) {
                            scales[i*parameters_per_site + h + searching_time + searching_m] = bottle_necks[j][k][1];
                        }else if(parameter_types[h] == selection_coeff){
                            scales[i*parameters_per_site + h + searching_time + searching_m] = bottle_necks[j][k][2];
                        }else if(parameter_types[h] == dominance_coeff){
                            scales[i*parameters_per_site + h + searching_time + searching_m] = bottle_necks[j][k][2] * 5;
                        }
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

                while((opt.max_value - opt.min_value) > bottle_necks[j][k][3] && opt.repeated_shrinkages < 4){
                    opt.iterate(to_be_optimized_function);

                    if(context.options.verbose_stderr){
                        cerr << "nelder-mead simplex_size: " << opt.simplex_size() << " output range: " << opt.max_value - opt.min_value << "\n";
                    }
                }

                



                found_parameters = opt.points[opt.max_index];


               // cout << "\n Result of search " << j + 1 << "/" << bottle_necks.size() << " " << k + 1 << "/" << bottle_necks[j].size() << " " << l + 1 << "/" << bottle_necks[j][k][0] << "\n";
                //cout << setprecision(15) << opt.max_value - context.fast_neutral_lnl<< "\n";
                //for(uint j = 0; j < found_parameters.size(); j++){
                    //cout << found_parameters[j] << "\n";
                //}
                //cout << "\n";

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










vector<double> search_sites(
    double chrom_size, 
    nelder_mead &opt,
    vector<vector<double>> sites,
    double width,
    double height,
    double depth
)
{
    vector<double> center_point(sites.size()*3);
    vector<double> scales(sites.size()*3);

    for(uint i = 0; i < sites.size(); i++) {
        center_point[i*3+0] = sites[i][0];
        center_point[i*3+1] = sites[i][1];
        center_point[i*3+2] = sites[i][2];
        scales[i*3+0] = width;
        scales[i*3+1] = height;
        scales[i*3+2] = height;
    }

    opt.populate_points(sites.size()*3, 1, center_point, scales);

    //random reflections
    for(uint i = 0; i < center_point.size(); i++){
        if(rand() < RAND_MAX/2){
            for(uint j = 0; j < opt.points.size(); j++){
                opt.points[j][i] = 2*center_point[i] - opt.points[j][i];
            }
        }
    }

    opt.enforce_bounds();
    
    opt.calculate_points(&to_be_optimized);

    while((opt.max_value - opt.min_value) > depth && opt.repeated_shrinkages < 4){
        
        opt.iterate(&to_be_optimized);
        
        if(context.options.verbose_stderr){
            cerr << "nelder-mead simplex_size: " << opt.simplex_size() << " output range: " << opt.max_value - opt.min_value << "\n";
        }
    }

    return opt.points[opt.max_index];
}



vector<double> search_sites_fast(double chrom_size, nelder_mead &opt, vector<vector<double>> sites, double width, double height, double depth){
    
    vector<double> center_point(sites.size()*3);
    vector<double> scales(sites.size()*3);

    for(uint i = 0; i < sites.size(); i++){
        center_point[i*3+0] = sites[i][0];
        center_point[i*3+1] = sites[i][1];
        center_point[i*3+2] = sites[i][2];
        scales[i*3+0] = width;
        scales[i*3+1] = height;
        scales[i*3+2] = height;
    }

    opt.populate_points(sites.size()*3, 1, center_point, scales);

    //random reflections
    for(uint i = 0; i < center_point.size(); i++){
        if(rand() < RAND_MAX/2){
            for(uint j = 0; j < opt.points.size(); j++){
                opt.points[j][i] = 2*center_point[i] - opt.points[j][i];
            }
        }
    }

    opt.enforce_bounds();

    opt.calculate_points(&to_be_optimized_only_near_sites);

    while((opt.max_value - opt.min_value) > depth && opt.repeated_shrinkages < 4){
        opt.iterate(&to_be_optimized_only_near_sites);
        
        if(context.options.verbose_stderr){
            cerr << "nelder-mead simplex_size: " << opt.simplex_size() << " output range: " << opt.max_value - opt.min_value << "\n";
        }
    }

    return opt.points[opt.max_index];
}


vector<double> search_sites_fast_fix_all_but_last(double chrom_size, nelder_mead &opt, vector<vector<double>> sites, double width, double height, double depth){
    
    vector<vector<double>> new_sites;

    vector<bool> in_new_sites(sites.size());

    for(uint i = 0; i < sites.size(); i++){

        in_new_sites[i] = false;

        double distance_from_last = abs(sites[i][0] - sites[sites.size() - 1][0]);

        if (distance_from_last <= fast_transitions_radius_in_morgans) {
            new_sites.push_back(sites[i]);
            in_new_sites[i] = true;
        }

    }
    

    vector<double> center_point(new_sites.size()*3);
    vector<double> scales(new_sites.size()*3);

    for(uint i = 0; i < new_sites.size(); i++){
        center_point[i*3+0] = new_sites[i][0];
        center_point[i*3+1] = new_sites[i][1];
        center_point[i*3+2] = new_sites[i][2];
        scales[i*3+0] = width;
        scales[i*3+1] = height;
        scales[i*3+2] = height;
    }
    

    opt.populate_points(new_sites.size()*3, 1, center_point, scales);

    //random reflections
    for(uint i = 0; i < center_point.size(); i++){
        if(rand() < RAND_MAX/2){
            for(uint j = 0; j < opt.points.size(); j++){
                opt.points[j][i] = 2*center_point[i] - opt.points[j][i];
            }
        }
    }
    
    
    opt.enforce_bounds();

    opt.calculate_points(&to_be_optimized_only_near_sites);
    

    while((opt.max_value - opt.min_value) > depth && opt.repeated_shrinkages < 4){
        
        opt.iterate(&to_be_optimized_only_near_sites);

        if(context.options.verbose_stderr){
            cerr << "nelder-mead simplex_size: " << opt.simplex_size() << " output range: " << opt.max_value - opt.min_value << "\n";
        }
    }
    
    vector<double> ans(sites.size()*3);
    
    int counter = 0;

    for(uint i = 0; i < sites.size(); i++){
        if(in_new_sites[i]){
            ans[i*3 + 0] = opt.points[opt.max_index][counter*3 + 0];
            ans[i*3 + 1] = opt.points[opt.max_index][counter*3 + 1];
            ans[i*3 + 2] = opt.points[opt.max_index][counter*3 + 2];

            counter++;
        }else{ 
            ans[i*3 + 0] = sites[i][0];
            ans[i*3 + 1] = sites[i][1];
            ans[i*3 + 2] = sites[i][2];
        }
    }
    
    return ans;
}







vector<double> search_sites_fast_dom0(double chrom_size, nelder_mead &opt, vector<vector<double>> sites, double width, double height, double depth){

    vector<double> center_point(sites.size()*2);
    vector<double> scales(sites.size()*2);

    for(uint i = 0; i < sites.size(); i++){
        center_point[i*2+0] = sites[i][0];
        center_point[i*2+1] = sites[i][1];
        scales[i*2+0] = width;
        scales[i*2+1] = height;
    }

    opt.populate_points(sites.size()*2, 1, center_point, scales);

    //random reflections
    for(uint i = 0; i < center_point.size(); i++){
        if(rand() < RAND_MAX/2){
            for(uint j = 0; j < opt.points.size(); j++){
                opt.points[j][i] = 2*center_point[i] - opt.points[j][i];
            }
        }
    }

    opt.enforce_bounds();

    opt.calculate_points(&to_be_optimized_pop0_dominant_fast);

    while((opt.max_value - opt.min_value) > depth && opt.repeated_shrinkages < 4){
        opt.iterate(&to_be_optimized_pop0_dominant_fast);

        if(context.options.verbose_stderr){
            cerr << "nelder-mead simplex_size: " << opt.simplex_size() << " output range: " << opt.max_value - opt.min_value << "\n";
        }
    }

    return opt.points[opt.max_index];
}

vector<double> search_sites_fast_dom1(double chrom_size, nelder_mead &opt, vector<vector<double>> sites, double width, double height, double depth){

    vector<double> center_point(sites.size()*2);
    vector<double> scales(sites.size()*2);

    for(uint i = 0; i < sites.size(); i++){
        center_point[i*2+0] = sites[i][0];
        center_point[i*2+1] = sites[i][1];
        scales[i*2+0] = width;
        scales[i*2+1] = height;
    }

    opt.populate_points(sites.size()*2, 1, center_point, scales);

    //random reflections
    for(uint i = 0; i < center_point.size(); i++){
        if(rand() < RAND_MAX/2){
            for(uint j = 0; j < opt.points.size(); j++){
                opt.points[j][i] = 2*center_point[i] - opt.points[j][i];
            }
        }
    }

    opt.enforce_bounds();

    opt.calculate_points(&to_be_optimized_pop1_dominant_fast);

    while((opt.max_value - opt.min_value) > depth && opt.repeated_shrinkages < 4){
        opt.iterate(&to_be_optimized_pop1_dominant_fast);

        if(context.options.verbose_stderr){
            cerr << "nelder-mead simplex_size: " << opt.simplex_size() << " output range: " << opt.max_value - opt.min_value << "\n";
        }
    }

    return opt.points[opt.max_index];
}


vector<double> search_sites_fast_additive(double chrom_size, nelder_mead &opt, vector<vector<double>> sites, double width, double height, double depth){

    vector<double> center_point(sites.size()*2);
    vector<double> scales(sites.size()*2);

    for(uint i = 0; i < sites.size(); i++){
        center_point[i*2+0] = sites[i][0];
        center_point[i*2+1] = sites[i][1];
        scales[i*2+0] = width;
        scales[i*2+1] = height;
    }

    opt.populate_points(sites.size()*2, 1, center_point, scales);

    //random reflections
    for(uint i = 0; i < center_point.size(); i++){
        if(rand() < RAND_MAX/2){
            for(uint j = 0; j < opt.points.size(); j++){
                opt.points[j][i] = 2*center_point[i] - opt.points[j][i];
            }
        }
    }

    opt.enforce_bounds();

    opt.calculate_points(&to_be_optimized_additive_fast);

    while((opt.max_value - opt.min_value) > depth && opt.repeated_shrinkages < 4){
        opt.iterate(&to_be_optimized_additive_fast);

        if(context.options.verbose_stderr){
            cerr << "nelder-mead simplex_size: " << opt.simplex_size() << " output range: " << opt.max_value - opt.min_value << "\n";
        }
    }

    return opt.points[opt.max_index];
}



vector<double> search_sites_fast_restricted_additive(double chrom_size, nelder_mead &opt, vector<vector<double>> sites, double width, double height, double depth){

    vector<double> center_point(sites.size());
    vector<double> scales(sites.size());

    for(uint i = 0; i < sites.size(); i++){
        center_point[i] = sites[i][0];
        scales[i] = height;
    }

    opt.populate_points(sites.size(), 1, center_point, scales);

    //random reflections
    for(uint i = 0; i < center_point.size(); i++){
        if(rand() < RAND_MAX/2){
            for(uint j = 0; j < opt.points.size(); j++){
                opt.points[j][i] = 2*center_point[i] - opt.points[j][i];
            }
        }
    }

    opt.enforce_bounds();

    opt.calculate_points(&to_be_optimized_restricted_additive_fast);

    while((opt.max_value - opt.min_value) > depth && opt.repeated_shrinkages < 4){
        opt.iterate(&to_be_optimized_restricted_additive_fast);

        if(context.options.verbose_stderr){
            cerr << "nelder-mead simplex_size: " << opt.simplex_size() << " output range: " << opt.max_value - opt.min_value << "\n";
        }
    }

    return opt.points[opt.max_index];
}




#endif
