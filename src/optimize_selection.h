#ifndef OPTIMIZE_SELECTION
#define OPTIMIZE_SELECTION
#include <vector>
#include <math.h>



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


    vector<double> smoothed_data_ancestry(data_ancestry.size());
    
    for(int i = 1; i < smoothed_data_ancestry.size() - 1; i++){

        double total = 0;
        double count = 0;
        
        for(int j = i; context.morgan_position[i] - context.morgan_position[j] < 0.001 && j >= 0; j--){
            total += data_ancestry[j];
            count ++;
        }
        for(int j = i + 1; context.morgan_position[j] - context.morgan_position[i] < 0.001 && j < smoothed_data_ancestry.size(); j++){
            total += data_ancestry[j];
            count ++;
        }
        
        smoothed_data_ancestry[i] = total/count;
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
            for (uint i = 0; i < selected_sites_count; i++) {
                cerr << "selection site: " << parameters[3*i + 0] << " with fitness: " << parameters[3*i + 1] << ",1," << parameters[3*i + 2] << "\n";
            }
        }else{
            for (uint i = 0; i < selected_sites_count; i++) {
                cerr << parameters[3*i + 0] << "\t" << parameters[3*i + 1] << "\t" << parameters[3*i + 2] << "\n";
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
            ss.site = parameters[i*3];
            ss.fitness.resize(3);
            ss.fitness[0] = parameters[i*3 + 1];
            ss.fitness[1] = 1;
            ss.fitness[2] = parameters[i*3 + 2];

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
    if(context.options.use_model_file)
        cores = 1;

    vector<mat> transition_matrices = calculate_transition_rates(
        context.n_recombs,
        selection_recomb_rates,
        fitnesses,
        context.options.m,
        context.options.generations,
        cores
    );
    
    if(parameters.size() == 0){
        last_calculated_transition_matricies = transition_matrices;
    }
    
    
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
    

    if(context.options.use_model_file)
        cores = 1;
    
    vector<mat> transition_matrices = alternative_fast_transition_rates (
        context.n_recombs,
        selection_recomb_rates,
        fitnesses,
        context.options.m,
        context.options.generations,
        cores
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














//each site is two parameters,p[0] p[1], which translates to (p[0], 1, p[1])
double to_be_optimized_pop0_dominant(vector<double> parameters){

    vector<double> new_parameters;

    int selected_sites_count = parameters.size() / 2;
    for(int i = 0; i < parameters.size(); i++){
        new_parameters.push_back(parameters[i]);
        if(i % 2 == 0){
            new_parameters.push_back(1);
        }
    }

    return to_be_optimized(new_parameters);
}

//each site is two parameters,p[0] p[1], which translates to (p[0], p[1], 1)
double to_be_optimized_pop1_dominant(vector<double> parameters){

    vector<double> new_parameters;

    int selected_sites_count = parameters.size() / 2;
    for(int i = 0; i < parameters.size(); i++){
        new_parameters.push_back(parameters[i]);
        if(i % 2 == 1){
            new_parameters.push_back(1);
        }
    }

    return to_be_optimized(new_parameters);
}

//each site is two parameters,p[0] p[1], which translates to (p[0], (1-p[1])/(1-p[1]/2), 1/(1-p[1]/2))
double to_be_optimized_additive(vector<double> parameters){

    vector<double> new_parameters;

    int selected_sites_count = parameters.size() / 2;
    for(int i = 0; i < parameters.size(); i++){
        if(i % 2 == 0){
            new_parameters.push_back(parameters[i]);
        }else{
            new_parameters.push_back((1-parameters[i])/(1-parameters[i]/2));
            new_parameters.push_back(1/(1-parameters[i]/2));
        }
        
    }

    return to_be_optimized(new_parameters);
}




//each site is two parameters,p[0] p[1], which translates to (p[0], 1, p[1])
double to_be_optimized_pop0_dominant_fast(vector<double> parameters){

    vector<double> new_parameters;

    int selected_sites_count = parameters.size() / 2;
    for(int i = 0; i < parameters.size(); i++){
        new_parameters.push_back(parameters[i]);
        if(i % 2 == 0){
            new_parameters.push_back(1);
        }
    }

    return to_be_optimized_only_near_sites(new_parameters);
}

//each site is two parameters,p[0] p[1], which translates to (p[0], p[1], 1)
double to_be_optimized_pop1_dominant_fast(vector<double> parameters){

    vector<double> new_parameters;

    int selected_sites_count = parameters.size() / 2;
    for(int i = 0; i < parameters.size(); i++){
        new_parameters.push_back(parameters[i]);
        if(i % 2 == 1){
            new_parameters.push_back(1);
        }
    }

    return to_be_optimized_only_near_sites(new_parameters);
}


//each site is two parameters,p[0] p[1], which translates to (p[0], (1-p[1])/(1-p[1]/2), 1/(1-p[1]/2))
double to_be_optimized_additive_fast(vector<double> parameters) {

    vector<double> new_parameters;

    int selected_sites_count = parameters.size() / 2;
    for(int i = 0; i < parameters.size(); i++){
        if(i % 2 == 0){
            new_parameters.push_back(parameters[i]);
        }else{
            new_parameters.push_back((1-parameters[i])/(1-parameters[i]/2));
            new_parameters.push_back(1/(1-parameters[i]/2));
        }
        
    }

    return to_be_optimized_only_near_sites(new_parameters);
}







double to_be_optimized_restricted(vector<double> parameters){

    vector<double> new_parameters;

    for(int i = 0; i < parameters.size()/2; i++){
        
        new_parameters.push_back(context.restricted_search_sites[i]);
        new_parameters.push_back(parameters[i*2 + 0]);
        new_parameters.push_back(parameters[i*2 + 1]);
        
    }

    return to_be_optimized(new_parameters);
}


double to_be_optimized_restricted_fast(vector<double> parameters) {

    vector<double> new_parameters;

    for(int i = 0; i < parameters.size()/2; i++){
        
        new_parameters.push_back(context.restricted_search_sites[i]);
        new_parameters.push_back(parameters[i*2 + 0]);
        new_parameters.push_back(parameters[i*2 + 1]);
        
    }

    return to_be_optimized_only_near_sites(new_parameters);
}












double to_be_optimized_restricted_additive(vector<double> parameters){

    vector<double> new_parameters;

    for(int i = 0; i < parameters.size(); i++){
        
        new_parameters.push_back(context.restricted_search_sites[i]);
        new_parameters.push_back((1-parameters[i])/(1-parameters[i]/2));
        new_parameters.push_back(1/(1-parameters[i]/2));
        
    }

    return to_be_optimized(new_parameters);
}

//each site is one parameters,p[0] which translates to (NA, (1-p[0])/(1-p[0]/2), 1/(1-p[0]/2))
double to_be_optimized_restricted_additive_fast(vector<double> parameters) {

    vector<double> new_parameters;

    for(int i = 0; i < parameters.size(); i++){
        
        new_parameters.push_back(context.restricted_search_sites[i]);
        new_parameters.push_back((1-parameters[i])/(1-parameters[i]/2));
        new_parameters.push_back(1/(1-parameters[i]/2));
        
    }

    return to_be_optimized_only_near_sites(new_parameters);
}







// (p[0]) => (NA, 1, p[0])
double to_be_optimized_restricted_dom0(vector<double> parameters){

    vector<double> new_parameters;

    for(int i = 0; i < parameters.size(); i++){
        
        new_parameters.push_back(context.restricted_search_sites[i]);
        new_parameters.push_back(1);
        new_parameters.push_back(parameters[i]);
        
    }

    return to_be_optimized(new_parameters);
}


double to_be_optimized_restricted_dom0_fast(vector<double> parameters) {

    vector<double> new_parameters;

    for(int i = 0; i < parameters.size(); i++){
        
        new_parameters.push_back(context.restricted_search_sites[i]);
        new_parameters.push_back(1);
        new_parameters.push_back(parameters[i]);
        
    }

    return to_be_optimized_only_near_sites(new_parameters);
}





// (p[0]) => (NA, p[0], 1)
double to_be_optimized_restricted_dom1(vector<double> parameters){

    vector<double> new_parameters;

    for(int i = 0; i < parameters.size(); i++){
        
        new_parameters.push_back(context.restricted_search_sites[i]);
        new_parameters.push_back(parameters[i]);
        new_parameters.push_back(1);
        
    }

    return to_be_optimized(new_parameters);
}


double to_be_optimized_restricted_dom1_fast(vector<double> parameters) {

    vector<double> new_parameters;

    for(int i = 0; i < parameters.size(); i++){
        
        new_parameters.push_back(context.restricted_search_sites[i]);
        new_parameters.push_back(parameters[i]);
        new_parameters.push_back(1);
        
    }

    return to_be_optimized_only_near_sites(new_parameters);
}
































typedef double (*lnl_function)(vector<double> parameters);


lnl_function to_be_optimized_variations (bool fast, bool restricted_site, bool additive, bool dom0, bool dom1) {

    

    if(restricted_site){
        if(fast){
            if(additive){
                return &to_be_optimized_restricted_additive_fast;
            }else if(dom0){
                return &to_be_optimized_restricted_dom0_fast;
            }else if(dom1){
                return &to_be_optimized_restricted_dom1_fast;
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
            }else{
                return &to_be_optimized;
            }
        }
    }
}





























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
            site[i] = parameters[j*parameters_per_site + i];
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
    vector<vector<double>> &sites,
    vector<vector<vector<double>>> &bottle_necks,
    double (*to_be_optimized_function) (vector<double>),
    vector<bool> is_loci,
    int parameters_per_site = 3
){
    vector<double> best_parameters;
    double best_ratio = -DBL_MAX;

    vector<double> found_parameters;
    

    for(int j = 0; j < bottle_necks.size(); j++){
        for(int k = 0; k < bottle_necks[j].size(); k++){
            for(int l = 0; l < bottle_necks[j][k][0]; l++){
                if(context.options.verbose_stderr){
                    cerr << "\n SEARCH " << j + 1 << "/" << bottle_necks.size() << " " << k + 1 << "/" << bottle_necks[j].size() << " " << l + 1 << "/" << bottle_necks[j][k][0] << "\n";
                }
                

                vector<double> center_point(sites.size()*parameters_per_site);
                vector<double> scales(sites.size()*parameters_per_site);

                
                for(uint i = 0; i < sites.size(); i++){
                    for(int h = 0; h < parameters_per_site; h++){
                        center_point[i*parameters_per_site + h] = sites[i][h];
                        if(is_loci[h]){
                            scales[i*parameters_per_site + h] = bottle_necks[j][k][1];
                        }else{
                            scales[i*parameters_per_site + h] = bottle_necks[j][k][2];
                        }
                    }
                }

                opt.populate_points(sites.size()*parameters_per_site, 1, center_point, scales);

                

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
                for(uint j = 0; j < found_parameters.size(); j++){
                    //cout << found_parameters[j] << "\n";
                }
                //cout << "\n";

                if(opt.max_value > best_ratio) {
                    best_ratio = opt.max_value;
                    best_parameters = found_parameters;
                }
            }
        }
        sites = parameters_to_sites(best_parameters, parameters_per_site);
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