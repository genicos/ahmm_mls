#ifndef SITE_EXAMINING
#define SITE_EXAMINING
#include <vector>
#include <math.h>

#include "optimize_selection.h"




vector<double> selection_opt::examine_sites(){
    vector<double> ans;


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







    
    double chrom_size = 0;
    for(uint i = 0; i < n_recombs.size(); i++){
        chrom_size += n_recombs[i];
    }

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

    vector<double> empty(0);

    context.neutral_lnl      = to_be_optimized(empty);
    context.fast_neutral_lnl = to_be_optimized_only_near_sites(empty);
    
    cerr << "\nNeutral likelihood: " << setprecision(15) << context.neutral_lnl << "\n";
    cerr << "\nFast neutral likelihood: " << setprecision(15) << context.fast_neutral_lnl << "\n";
    cout << "\nNeutral likelihood: " << setprecision(15) << context.neutral_lnl << "\n";
    cout << "\nFast neutral likelihood: " << setprecision(15) << context.fast_neutral_lnl << "\n";

    neutral_transition_matrices = last_calculated_transition_matricies;


    cerr << setprecision(15);
    cout << setprecision(15);



    for(int i = 0; i < options.site_file_morgan_positions.size(); i++) {
        
        string option = options.site_file_options[i][0];
        

        
        char restrict_site_op = 's';
        bool restrict_site = option.find(restrict_site_op) != string::npos;
        
        char additive_op = 'a';
        bool additive = option.find(additive_op) != string::npos;

        char dom0_op = 'D';
        bool dom0 = option.find(dom0_op) != string::npos;
        
        char dom1_op = 'd';
        bool dom1 = option.find(dom1_op) != string::npos;

        char sel0_op = '0';
        bool sel0 = option.find(sel0_op) != string::npos;

        char sel1_op = '1';
        bool sel1 = option.find(sel1_op) != string::npos;


        char nsch_op = 'n';
        bool nsch = option.find(nsch_op) != string::npos;



        int site_count = options.site_file_morgan_positions[i].size();
        int parameters_per_site = ( !restrict_site + !(additive || dom0 || dom1) + 1);
        int parameter_count = parameters_per_site * site_count;


        vector<bool> is_loci(0);
        if(!restrict_site){
            is_loci.push_back(true);
        }else{
            context.restricted_search_sites.resize(site_count);
            for(int j = 0; j < site_count; j++) {
                context.restricted_search_sites[j] = options.site_file_morgan_positions[i][j];
            }
            
        }

        is_loci.push_back(false);
        if(!(additive || dom0 || dom1)){
            is_loci.push_back(false);
        }



        if(nsch) { //If we have the no search option, then we skip to directly computing the likelihood
            vector<double> parameters(0);

            for(int j = 0; j < site_count; j++){
                if(!restrict_site){ 
                    parameters.push_back(options.site_file_morgan_positions[i][j]);
                }

                if(additive || dom0 || dom1){
                    parameters.push_back(options.site_file_low_bounds[i][j]);
                }else{
                    parameters.push_back(options.site_file_low_bounds[i][j]);
                    parameters.push_back(options.site_file_high_bounds[i][j]);
                }
            }

            vector<vector<double>> sites = parameters_to_sites(parameters, parameters_per_site);

            for(int j = 0; j < site_count; j++){
                
                int index = 0;

                cout << "site:\t";
                cerr << "site:\t";

                if(restrict_site){
                    cout << context.restricted_search_sites[j] << "\t";
                    cerr << context.restricted_search_sites[j] << "\t";
                }else{
                    cout << sites[j][index] << "\t";
                    cerr << sites[j][index] << "\t";
                    index ++;
                }

                if(additive){
                    cout << (1-sites[j][index]) << ","<< (1-sites[j][index]/2) << ",1\n";
                    cerr << (1-sites[j][index]) << ","<< (1-sites[j][index]/2) << ",1\n";
                }else if(dom0){
                    cout << "1,1," << sites[j][index] << "\n";
                    cerr << "1,1," << sites[j][index] << "\n";
                }else if(dom1){
                    cout << sites[j][index] << ",1,1\n";
                    cerr << sites[j][index] << ",1,1\n";
                }else{
                    cout << sites[j][index] << ",1," << sites[j][index + 1] << "\n";
                    cerr << sites[j][index] << ",1," << sites[j][index + 1] << "\n";
                }
            }

            double fast_singular_lnl = to_be_optimized_variations(true, restrict_site, additive, dom0, dom1) (parameters);

            cout << setprecision(15) << "fast lnl ratio:\t" << fast_singular_lnl - context.fast_neutral_lnl<< "\n";
            cerr << setprecision(15) << "fast lnl ratio:\t" << fast_singular_lnl - context.fast_neutral_lnl<< "\n";

            continue;
        }
        

        vector<vector<double>> sites;
        
        vector<double> singular_starting_parameters(0);

        for(int j = 0; j < site_count; j++){
            if(!restrict_site){ 
                singular_starting_parameters.push_back(options.site_file_morgan_positions[i][j]);
            }

            if(additive){
                singular_starting_parameters.push_back(0);
            }else if (dom0 || dom1){
                singular_starting_parameters.push_back(1);
            }else{
                singular_starting_parameters.push_back(1);
                singular_starting_parameters.push_back(1);
            }
        }
        

        // single site additive selection testing
        sites = parameters_to_sites(singular_starting_parameters, parameters_per_site);

        

        optimizer.init_bounds(parameter_count, 0.001);


        for(int j = 0; j < site_count; j++){
            if(!restrict_site){
                optimizer.min_bounds[j*parameters_per_site] = options.site_file_low_bounds[i][j];
                optimizer.max_bounds[j*parameters_per_site] = options.site_file_high_bounds[i][j];
            }

            if(additive){
                if(sel0){
                    optimizer.max_bounds[j*parameters_per_site + (!restrict_site)] = 0;
                }else if(sel1){
                    optimizer.min_bounds[j*parameters_per_site + (!restrict_site)] = 0;
                    optimizer.max_bounds[j*parameters_per_site + (!restrict_site)] = 1;
                }else{
                    optimizer.max_bounds[j*parameters_per_site + (!restrict_site)] = 1;
                }
            }else if (dom0){
                if(sel0){
                    optimizer.min_bounds[j*parameters_per_site + (!restrict_site)] = 0;
                    optimizer.max_bounds[j*parameters_per_site + (!restrict_site)] = 1;
                }else if(sel1){
                    optimizer.min_bounds[j*parameters_per_site + (!restrict_site)] = 1;
                }else{
                    optimizer.min_bounds[j*parameters_per_site + (!restrict_site)] = 0;
                }
            }else if (dom1){
                if(sel0){
                    optimizer.min_bounds[j*parameters_per_site + (!restrict_site)] = 1;
                }else if(sel1){
                    optimizer.min_bounds[j*parameters_per_site + (!restrict_site)] = 0;
                    optimizer.max_bounds[j*parameters_per_site + (!restrict_site)] = 1;
                }else{
                    optimizer.min_bounds[j*parameters_per_site + (!restrict_site)] = 0;
                }
            }else{
                optimizer.min_bounds[j*parameters_per_site + (!restrict_site)] = 0;
                optimizer.min_bounds[j*parameters_per_site + (!restrict_site) + 1] = 0;
            }
        }


        vector<double> singular_optimized_parameters = multi_level_optimization(
            chrom_size,
            optimizer,
            sites,
            bottle_necks,
            to_be_optimized_variations(true, restrict_site, additive, dom0, dom1),
            is_loci,
            parameters_per_site
        );


        for(int j = 0; j < site_count; j++){
            
            int index = 0;

            cout << "site:\t";
            cerr << "site:\t";

            if(restrict_site){
                cout << context.restricted_search_sites[j] << "\t";
                cerr << context.restricted_search_sites[j] << "\t";
            }else{
                cout << sites[j][index] << "\t";
                cerr << sites[j][index] << "\t";
                index ++;
            }

            if(additive){
                cout << (1-sites[j][index]) << ","<< (1-sites[j][index]/2) << ",1\n";
                cerr << (1-sites[j][index]) << ","<< (1-sites[j][index]/2) << ",1\n";
            }else if(dom0){
                cout << "1,1," << sites[j][index] << "\n";
                cerr << "1,1," << sites[j][index] << "\n";
            }else if(dom1){
                cout << sites[j][index] << ",1,1\n";
                cerr << sites[j][index] << ",1,1\n";
            }else{
                cout << sites[j][index] << ",1," << sites[j][index + 1] << "\n";
                cerr << sites[j][index] << ",1," << sites[j][index + 1] << "\n";
            }
        }
        

        double fast_singular_lnl = to_be_optimized_variations(true, restrict_site, additive, dom0, dom1) (singular_optimized_parameters);
        
        cout << setprecision(15) << "fast lnl ratio:\t" << fast_singular_lnl - context.fast_neutral_lnl<< "\n";
        cerr << setprecision(15) << "fast lnl ratio:\t" << fast_singular_lnl - context.fast_neutral_lnl<< "\n";



        //TODO command line option to turn this off!!!!
        /*
        double singular_lnl = to_be_optimized_variations(false, restrict_site, additive, dom0, dom1) (singular_optimized_parameters);

        cout << setprecision(15) << "lnl ratio:\t" << singular_lnl - context.neutral_lnl << "\n";
        cerr << setprecision(15) << "lnl ratio:\t" << singular_lnl - context.neutral_lnl << "\n";
        */
        
    }

    return ans;
}

#endif
