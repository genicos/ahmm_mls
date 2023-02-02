#ifndef __EXAMINE_MODELS_H
#define __EXAMINE_MODELS_H
#include <vector>
#include <math.h>

#include "optimize_selection.h"

void selection_opt::examine_models() {

    // Defining multi-level optimization parameters //
    vector<vector<vector<double>>> bottle_necks;    //

    vector<vector<double>> shallow;
    
    vector<double> shallow_short(4);
    shallow_short[0] = 5;
    shallow_short[1] = 0.03;
    shallow_short[2] = 0.01;
    shallow_short[3] = 5;

    vector<double> shallow_tall(4);
    shallow_tall[0] = 5;
    shallow_tall[1] = 0.03;
    shallow_tall[2] = 0.05;
    shallow_tall[3] = 5;

    shallow.push_back(shallow_short);
    shallow.push_back(shallow_tall);

    vector<vector<double>> deep;

    vector<double> deep_short(4);
    deep_short[0] = 5;
    deep_short[1] = 0.01;
    deep_short[2] = 0.005;
    deep_short[3] = 1;

    vector<double> deep_tall(4);
    deep_tall[0] = 5;
    deep_tall[1] = 0.01;
    deep_tall[2] = 0.01;
    deep_tall[3] = 1;

    deep.push_back(deep_short);
    deep.push_back(deep_tall);

    bottle_necks.push_back(shallow);
    bottle_necks.push_back(deep);                   //
    //////////////////////////////////////////////////


    double nelder_mead_reflection  = 1;
    double nelder_mead_contraction = 0.61803398875;
    double nelder_mead_expansion   = 2;
    double nelder_mead_shrinkage   = 0.5;

    nelder_mead optimizer(
        nelder_mead_reflection,
        nelder_mead_contraction,
        nelder_mead_expansion,
        nelder_mead_shrinkage
    );






    context = *this;

    /*
    vector<double> empty(0);
    context.neutral_lnl      = to_be_optimized(empty);
    context.fast_neutral_lnl = to_be_optimized_only_near_sites(empty);
    
    cerr << "\nNeutral likelihood: " << setprecision(15) << context.neutral_lnl << "\n";
    cerr << "\nFast neutral likelihood: " << setprecision(15) << context.fast_neutral_lnl << "\n";
    cout << "Neutral likelihood: " << setprecision(15) << context.neutral_lnl << "\n";
    cout << "Fast neutral likelihood: " << setprecision(15) << context.fast_neutral_lnl << "\n";
    neutral_transition_matrices = last_calculated_transition_matricies;
    */


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

            double lnl = general_to_be_optimized_fast(initial_parameters);

            cout << setprecision(15) << "fast lnl ratio:\t" << lnl - context.fast_neutral_lnl << "\n";
            cerr << setprecision(15) << "fast lnl ratio:\t" << lnl - context.fast_neutral_lnl << "\n";
            
        
        
        
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
                    if(options.mls_searches[i].min_bound_l[j] == -1)
                        optimizer.min_bounds[p] = 0;
                    else
                        optimizer.min_bounds[p] = options.mls_searches[i].min_bound_l[j];
                    
                    if(options.mls_searches[i].max_bound_l[j] == -1)
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
        

            vector<double> optimized_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                initial_parameters,
                bottle_necks,
                &general_to_be_optimized_fast,
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

            


            double lnl = general_to_be_optimized_fast (optimized_parameters);
            
            cout << setprecision(15) << "fast lnl ratio:\t" << lnl - context.fast_neutral_lnl << "\n";
            cerr << setprecision(15) << "fast lnl ratio:\t" << lnl - context.fast_neutral_lnl << "\n";


            initial_parameters = optimized_parameters;


        }


        char model_posterior_op = 'm';
        bool model_posterior_printing = options.mls_searches[i].posterior_options.find(model_posterior_op) != string::npos;

        char data_posterior_op = 'd';
        bool data_posterior_printing = options.mls_searches[i].posterior_options.find(data_posterior_op) != string::npos;


        if (model_posterior_printing || data_posterior_printing)
            general_to_be_optimized (initial_parameters);

        if ( model_posterior_printing ) {

            ofstream model_ancestry;

            string file_name = "model_ancestry_" + options.mls_searches[i].search_name + ".tsv";

            model_ancestry.open(file_name);

            model_ancestry << "morgan_pos\tlocal_ancestry\n";

            for(int i = 1; i < n_recombs.size() - 1; i++){

                model_ancestry << morgan_position[i]  << "\t" << local_ancestries[i] << "\n";
            }
                
            model_ancestry.close();
        }

    }
}

#endif