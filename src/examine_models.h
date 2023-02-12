#ifndef __EXAMINE_MODELS_H
#define __EXAMINE_MODELS_H
#include <vector>
#include <math.h>

#include "optimize_selection.h"

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
            
            sample_ancestry[i] = conv_to< vector<double> >::from(smoothed_probs);
            
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




void selection_opt::examine_models() {

    // Defining optimization parameters
    vector<vector<vector<double>>> search_restarts;   

    vector<vector<double>> shallow;
    
    vector<double> shallow_short = {5, 0.03, 0.01, 5};
    vector<double> shallow_tall  = {5, 0.03, 0.05, 5};
    shallow.push_back(shallow_short);
    shallow.push_back(shallow_tall);

    vector<vector<double>> deep;

    vector<double> deep_short = {5, 0.01, 0.005, 1};
    vector<double> deep_tall  = {5, 0.01, 0.01 , 1};
    deep.push_back(deep_short);
    deep.push_back(deep_tall);

    search_restarts.push_back(shallow);
    search_restarts.push_back(deep);

    double nelder_mead_reflection  = 1;
    double nelder_mead_contraction = 0.5;
    double nelder_mead_expansion   = 2;
    double nelder_mead_shrinkage   = 0.5;

    nelder_mead optimizer(
        nelder_mead_reflection,
        nelder_mead_contraction,
        nelder_mead_expansion,
        nelder_mead_shrinkage
    );


    context = *this;



    cerr << setprecision(15);
    cout << setprecision(15);

    for(int i = 0; i < options.mls_searches.size(); i++){

        cerr << options.mls_searches[i].search_string << "\n";
        cout << options.mls_searches[i].search_string << "\n";

        setup_searches(options.mls_searches[i]);

        
        vector<parameter_type> parameter_types(0);

        
        int site_count = options.mls_searches[i].search_l.size(); 
        int parameter_count = 0;

        if(options.mls_searches[i].search_m){
            parameter_types.push_back(admix_frac);
            parameter_count++;
        }
        if(options.mls_searches[i].search_t){
            parameter_types.push_back(time_since_admix);
            parameter_count++;
        }

        for(int j = 0; j < site_count; j++){
            if(options.mls_searches[i].search_l[j]){
                parameter_types.push_back(location);
                parameter_count++;
            }
            if(options.mls_searches[i].search_h[j]){
                parameter_types.push_back(dominance_coeff);
                parameter_count++;
            }
            if(options.mls_searches[i].search_s[j]){
                parameter_types.push_back(selection_coeff);
                parameter_count++;
            }
        }

        vector<double> initial_parameters(parameter_count);


        if(parameter_count == 0) { //No searching, just fitting a model

            double lnl = to_be_optimized_fast(initial_parameters);

            cout << setprecision(15) << "fast lnl:\t" << lnl << "\n";
            cerr << setprecision(15) << "fast lnl:\t" << lnl << "\n";
        
        } else { // Searching

            optimizer.init_bounds(parameter_count, 0.001);
            
            int p = 0;
            if (options.mls_searches[i].search_m) {
                optimizer.min_bounds[p] = 0;
                optimizer.max_bounds[p] = 1;
                initial_parameters[p] = options.mls_searches[i].start_m;
                p++;
            }
            
            if (options.mls_searches[i].search_t) {
                optimizer.min_bounds[p] = 1;
                optimizer.max_bounds[p] = 10000;
                initial_parameters[p] = options.mls_searches[i].start_t;
                p++;
            }


            for(int j = 0; j < site_count; j++) {
                
                if(options.mls_searches[i].search_l[j]) {
                    if(options.mls_searches[i].min_bound_l[j] < 0)
                        optimizer.min_bounds[p] = 0;
                    else
                        optimizer.min_bounds[p] = options.mls_searches[i].min_bound_l[j];
                    
                    if(options.mls_searches[i].max_bound_l[j] == -1 || options.mls_searches[i].max_bound_l[j] > chrom_size)
                        optimizer.max_bounds[p] = chrom_size;
                    else
                        optimizer.max_bounds[p] = options.mls_searches[i].max_bound_l[j];
                    
                    initial_parameters[p] = options.mls_searches[i].start_l[j];
                    p++;
                }

                if(options.mls_searches[i].search_h[j]) { //TODO should h always be [0,1]?
                    optimizer.min_bounds[p] = 0;
                    optimizer.max_bounds[p] = 1;
                    initial_parameters[p] = options.mls_searches[i].start_h[j];
                    p++;
                }

                if(options.mls_searches[i].search_s[j]) {
                    if( options.mls_searches[i].s_pos[j] )
                        optimizer.min_bounds[p] = 0;
                    
                    if( options.mls_searches[i].s_neg[j] )
                        optimizer.max_bounds[p] = 0;
                    else
                        optimizer.max_bounds[p] = 1;
                    
                    initial_parameters[p] = options.mls_searches[i].start_s[j];
                    p++;
                }
            }

            for(int h = 0; h < initial_parameters.size(); h++){
                cerr << initial_parameters[h] << "\n";
                cerr << optimizer.min_bounds[h] << "\n";
                cerr << optimizer.max_bounds[h] << "\n\n\n";
            }
        

            vector<double> optimized_parameters = multi_restart_optimization(
                chrom_size,
                optimizer,
                initial_parameters,
                search_restarts,
                &to_be_optimized_fast,
                parameter_types
            );

            cerr << "\n\nBest Model\n";

            

            vector<double> parameters = convert_parameters_to_long_form(optimized_parameters);

            if (global_search.search_m) {
                cerr << "m: " << parameters[0] << "\n";
                cout << "m: " << parameters[0] << "\n";
            }
            if (global_search.search_t) {
                cerr << "Time:    " << parameters[global_search.search_m] << "\n";
                cout << "Time:    " << parameters[global_search.search_m] << "\n";
            }
            for (uint i = 0; i < site_count; i++) {
                cerr << "site: " << parameters[3*i + 0 + global_search.search_m + global_search.search_t] << " with fitness: " << parameters[3*i + 1 + global_search.search_m + global_search.search_t] << ",1," << parameters[3*i + 2 + global_search.search_m + global_search.search_t] << "\n";
                cout << "site: " << parameters[3*i + 0 + global_search.search_m + global_search.search_t] << " with fitness: " << parameters[3*i + 1 + global_search.search_m + global_search.search_t] << ",1," << parameters[3*i + 2 + global_search.search_m + global_search.search_t] << "\n";
            }

            


            double lnl = to_be_optimized_fast (optimized_parameters);
            
            cout << setprecision(15) << "fast lnl ratio:\t" << lnl << "\n";
            cerr << setprecision(15) << "fast lnl ratio:\t" << lnl << "\n";


            initial_parameters = optimized_parameters;


        }


        char model_posterior_op = 'm';
        bool model_posterior_printing = options.mls_searches[i].posterior_options.find(model_posterior_op) != string::npos;

        char data_posterior_op = 'd';
        bool data_posterior_printing = options.mls_searches[i].posterior_options.find(data_posterior_op) != string::npos;

        char sample_posterior_op = 's';
        bool sample_posterior_printing = options.mls_searches[i].posterior_options.find(sample_posterior_op) != string::npos;


        if (model_posterior_printing || data_posterior_printing || sample_posterior_op)
            to_be_optimized (initial_parameters);



        if ( model_posterior_printing ) {

            ofstream model_ancestry;

            string file_name = "model_ancestry_" + options.mls_searches[i].search_name + ".tsv";

            model_ancestry.open(file_name);

            model_ancestry << "position\tposition_morg\tlocal_ancestry";

            for(int i = 1; i < n_recombs.size(); i++) {

                model_ancestry << "\n" << context.position[i] << "\t" << morgan_position[i]  << "\t" << local_ancestries[i];
            }
                
            model_ancestry.close();
        }



        if ( sample_posterior_printing ) {
            get_local_ancestry(last_calculated_transition_matricies);



            pulse first_pulse;
            first_pulse.time = 10000000; //TODO if im allowing time to be searched, change this
            first_pulse.time_fixed = true;
            first_pulse.type = 1;
            first_pulse.proportion = 1 - options.mls_searches[i].start_m;
            first_pulse.proportion_fixed = true;
            first_pulse.entry_order = 0;

            pulse second_pulse;
            second_pulse.time = options.mls_searches[i].start_t;
            second_pulse.time_fixed = true;
            second_pulse.type = 0;
            second_pulse.proportion = options.mls_searches[i].start_m;
            second_pulse.proportion_fixed = true;
            second_pulse.entry_order = 1;

            vector<pulse> optimum(2);
            optimum[0] = first_pulse;
            optimum[1] = second_pulse;

            
            cerr << "forward-backward posterior decoding and printing\t\t\t";
            for ( int m = 0 ; m < markov_chain_information.size() ; m ++ ) {
                context.markov_chain_information[m].combine_prob( context.position, context.state_list, context.chromosomes, options.output_pulses, optimum ) ;
            }
            
        }



        if ( data_posterior_printing ) {

            vector<double> data_local_ancestry = get_local_ancestry(last_calculated_transition_matricies);

            ofstream data_ancestry;

            string file_name = "data_ancestry_" + options.mls_searches[i].search_name + ".tsv";

            data_ancestry.open(file_name);

            data_ancestry << "position\tposition_morg\tlocal_ancestry";

            for(int i = 1; i < n_recombs.size(); i++) {

                data_ancestry << "\n" << context.position[i] << "\t" << morgan_position[i]  << "\t" << data_local_ancestry[i];
            }

            data_ancestry.close();

        }

    }
}

#endif