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
    shallow_short[3] = 20;

    vector<double> shallow_tall(4);
    shallow_tall[0] = 5;
    shallow_tall[1] = 0.03;
    shallow_tall[2] = 0.05;
    shallow_tall[3] = 20;

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


    // Defining multi-level optimization parameters ///////
    vector<vector<vector<double>>> two_site_bottle_necks;//

    
    shallow_short[0] = 5;
    shallow_short[1] = 0.005;
    shallow_short[2] = 0.005;
    shallow_short[3] = 5;

    
    shallow_tall[0] = 5;
    shallow_tall[1] = 0.005;
    shallow_tall[2] = 0.03;
    shallow_tall[3] = 5;

    shallow[0] = shallow_short;
    shallow[1] = shallow_tall;


    deep_short[0] = 5;
    deep_short[1] = 0.002;
    deep_short[2] = 0.002;
    deep_short[3] = 1;

    deep_tall[0] = 5;
    deep_tall[1] = 0.002;
    deep_tall[2] = 0.005;
    deep_tall[3] = 1;

    deep[0] = deep_short;
    deep[1] = deep_tall;

    two_site_bottle_necks.push_back(shallow);
    two_site_bottle_necks.push_back(deep);               //
    ///////////////////////////////////////////////////////





    
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

    context.neutral_lnl = to_be_optimized(empty);
    neutral_lnl = context.neutral_lnl;                  //TODO i dont know where this is declared?
    
    cerr << "\nNeutral likelihood: " << setprecision(15) << neutral_lnl << "\n";

    neutral_transition_matrices = last_calculated_transition_matricies;


    cerr << setprecision(15);
    cout << setprecision(15);

    for(int i = 0; i < options.site_file_morgan_positions.size(); i++){

        double bound_size = options.site_file_high_bounds[i] 
                - options.site_file_low_bounds[i];
        
        if (options.site_file_options[i].compare("o") == 0) {

            vector<double> starting_parameters(3);
            starting_parameters[0] = options.site_file_morgan_positions[i];
            starting_parameters[1] = 1;
            starting_parameters[2] = 1;

            starting_parameters = grid_search(starting_parameters, 0.001, 0.05, 3, 6);

            if(options.verbose_stderr){
                cerr << "Best grid search result:\n";
                cerr << starting_parameters[0] << "\n";
                cerr << starting_parameters[1] << "\n";
                cerr << starting_parameters[2] << "\n";
            }


            vector<vector<double>> sites = parameters_to_sites(starting_parameters);

            
            optimizer.init_bounds(3, min(0.01, bound_size/2) );
            optimizer.min_bounds[0] = options.site_file_low_bounds[i];
            optimizer.max_bounds[0] = options.site_file_high_bounds[i];

            optimizer.min_bounds[1] = 0;
            optimizer.min_bounds[2] = 0;

            
            vector<double> best_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast
            );

                
            double lnl = to_be_optimized(best_parameters);

            cout << setprecision(15) << "\n\nUnrestricted optimization lnl:\t" << lnl << "\n";
            cout << "\n\nneutral lnl:\t" << neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << sites[0][1] << ",1," << sites[0][2] << " or " << sites[0][1]/max(sites[0][1], sites[0][2]) << ","<< 1/max(sites[0][1], sites[0][2]) <<"," << sites[0][2]/max(sites[0][1], sites[0][2]) << "\n";
            
            if(options.verbose_stderr){
                cerr << "\n\nUnrestricted optimization lnl:\t" << lnl << "\n";
                cerr << "\n\nneutral lnl:\t" << neutral_lnl << "\n";
                cerr << "site:\t" << sites[0][0] << "\t" << sites[0][1] << ",1," << sites[0][2] << " or " << sites[0][1]/max(sites[0][1], sites[0][2]) << ","<< 1/max(sites[0][1], sites[0][2]) <<"," << sites[0][2]/max(sites[0][1], sites[0][2]) << "\n";
            }else{
                cerr << "FINAL SITE:\n";
                cerr << "Neutral lnl\t" << setprecision(15) << neutral_lnl << "\tNeutral lnl\t" << lnl << "\n";
                cerr << "site:\t" << sites[0][0] << "\t" << sites[0][1] << ",1," << sites[0][2] << " or " << sites[0][1]/max(sites[0][1], sites[0][2]) << ","<< 1/max(sites[0][1], sites[0][2]) <<"," << sites[0][2]/max(sites[0][1], sites[0][2]) << "\n";
            }
        }
        
        else if (options.site_file_options[i].compare("d") == 0) {

            //Dominance testing
            vector<double> dom_starting_parameters(2);
            dom_starting_parameters[0] = options.site_file_morgan_positions[i];
            dom_starting_parameters[1] = 1;

            vector<double> add_starting_parameters(2);
            add_starting_parameters[0] = options.site_file_morgan_positions[i];
            add_starting_parameters[1] = 0;
            


            // Pop0 dominant testing
            vector<double> dom0_starting_parameters = grid_search_dominant0(dom_starting_parameters, 0.001, 0.05, 3, 6);

            if(options.verbose_stderr){
                cerr << "Best grid search result:\n";
                cerr << dom0_starting_parameters[0] << "\n";
                cerr << dom0_starting_parameters[1] << "\n";
            }

            vector<vector<double>> sites = parameters_to_sites(dom0_starting_parameters, 2);
            

            optimizer.init_bounds(2, min(0.01, bound_size/2) );
            optimizer.min_bounds[0] = options.site_file_low_bounds[i];
            optimizer.max_bounds[0] = options.site_file_high_bounds[i];

            optimizer.min_bounds[1] = 0;
            
            vector<double> pop0_dom_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast_dom0,
                2
            );

            double dom0_lnl = to_be_optimized_pop0_dominant(pop0_dom_parameters);

            cout << setprecision(15) << "\n\nPop0 dom lnl ratio: " << dom0_lnl - neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << "1,1," << sites[0][1] << " or " << 1/max(1.0, sites[0][1]) << "," << 1/max(1.0, sites[0][1]) <<"," << sites[0][1]/max(1.0, sites[0][1]) << "\n";

            cerr << "\n\nPop0 dom lnl ratio: " << dom0_lnl - neutral_lnl << "\n";
            cerr << "site:\t" << sites[0][0] << "\t" << "1,1," << sites[0][1] << " or " << 1/max(1.0, sites[0][1]) << "," << 1/max(1.0, sites[0][1]) <<"," << sites[0][1]/max(1.0, sites[0][1]) << "\n";




            // Pop1 dominant testing
            vector<double> dom1_starting_parameters = grid_search_dominant1(dom_starting_parameters, 0.001, 0.05, 3, 6);

            if(options.verbose_stderr) {
                cerr << "Best grid search result:\n";
                cerr << dom1_starting_parameters[0] << "\n";
                cerr << dom1_starting_parameters[1] << "\n";
            }

            sites = parameters_to_sites(dom1_starting_parameters, 2);
            
            optimizer.init_bounds(2, min(0.01, bound_size/5) );
            optimizer.min_bounds[0] = options.site_file_low_bounds[i];
            optimizer.max_bounds[0] = options.site_file_high_bounds[i];

            optimizer.min_bounds[1] = 0;

            vector<double> pop1_dom_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast_dom1,
                2
            );

            double dom1_lnl = to_be_optimized_pop1_dominant(pop1_dom_parameters);

            cout << setprecision(15) << "\n\nPop1 dom lnl ratio:\t" << dom1_lnl - neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << sites[0][1] << ",1,1  or " << sites[0][1]/max(1.0, sites[0][1]) << "," << 1/max(1.0, sites[0][1]) <<"," << 1/max(1.0, sites[0][1]) << "\n";

            cerr << "\n\nPop1 dom lnl ratio:\t" << dom1_lnl - neutral_lnl << "\n";
            cerr << "site:\t" << sites[0][0] << "\t" << sites[0][1] << ",1,1  or " << sites[0][1]/max(1.0, sites[0][1]) << "," << 1/max(1.0, sites[0][1]) <<"," << 1/max(1.0, sites[0][1]) << "\n";



            // Additive selection testing
            add_starting_parameters = grid_search_additive(add_starting_parameters, 0.001, 0.05, 3, 6);

            if(options.verbose_stderr){
                cerr << "Best grid search result:\n";
                cerr << add_starting_parameters[0] << "\n";
                cerr << add_starting_parameters[1] << "\n";
            }

            sites = parameters_to_sites(add_starting_parameters, 2);

            optimizer.init_bounds(2, min(0.01, bound_size/5) );
            optimizer.min_bounds[0] = options.site_file_low_bounds[i];
            optimizer.max_bounds[0] = options.site_file_high_bounds[i];

            optimizer.min_bounds[1] = -1;

            vector<double> additive_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast_additive,
                2
            );

            double add_lnl = to_be_optimized_additive(additive_parameters);

            cout << setprecision(15) << "\n\nAdditive lnl ratio:\t" << add_lnl - neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << (1-sites[0][1]) << ","<< (1-sites[0][1]/2) << ",1\t1," << (1-sites[0][1]/2)/(1-sites[0][1]) << "," << 1/(1-sites[0][1]) << "\n";

            cerr << "\n\nAdditive lnl ratio:\t" << add_lnl - neutral_lnl << "\n";
            cerr << "site:\t" << sites[0][0] << "\t" << (1-sites[0][1]) << ","<< (1-sites[0][1]/2) << ",1\t1," << (1-sites[0][1]/2)/(1-sites[0][1]) << "," << 1/(1-sites[0][1]) << "\n";


            cerr << "\n\nCompleted Dominance testing of site at " << options.site_file_positions[i] << "\n";
            cerr << setprecision(15) << "Pop0 dominant lnl:\t" << dom0_lnl << "\n";
            cerr << setprecision(15) << "Pop1 dominant lnl:\t" << dom1_lnl << "\n";
            cerr << setprecision(15) << "Additive selection lnl:\t" << add_lnl << "\n";

            cout << "\n\nCompleted Dominance testing of site at " << options.site_file_positions[i] << "\n";
            cout << setprecision(15) << "Pop0 dominant lnl:\t" << dom0_lnl << "\n";
            cout << setprecision(15) << "Pop1 dominant lnl:\t" << dom1_lnl << "\n";
            cout << setprecision(15) << "Additive selection lnl:\t" << add_lnl << "\n";
        }

        else if (options.site_file_options[i].compare("T") == 0){

            vector<vector<double>> sites;


            
            vector<double> add_starting_parameters(2);
            add_starting_parameters[0] = options.site_file_morgan_positions[i];
            add_starting_parameters[1] = 0;

            // single site additive selection testing
            sites = parameters_to_sites(add_starting_parameters, 2);


            optimizer.init_bounds(2, min(0.01, bound_size/5) );
            optimizer.min_bounds[0] = options.site_file_low_bounds[i];
            optimizer.max_bounds[0] = options.site_file_high_bounds[i];

            optimizer.min_bounds[1] = -1;

            vector<double> additive_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast_additive,
                2
            );

            double add_lnl = to_be_optimized_additive(additive_parameters);


            cout << setprecision(15) << "\n\nAdditive lnl ratio:\t" << add_lnl - neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << (1-sites[0][1]) << ","<< (1-sites[0][1]/2) << ",1\t1," << (1-sites[0][1]/2)/(1-sites[0][1]) << "," << 1/(1-sites[0][1]) << "\n";

            cerr << setprecision(15) << "\n\nAdditive lnl ratio:\t" << add_lnl - neutral_lnl << "\n";
            cerr << "site:\t" << sites[0][0] << "\t" << (1-sites[0][1]) << ","<< (1-sites[0][1]/2) << ",1\t1," << (1-sites[0][1]/2)/(1-sites[0][1]) << "," << 1/(1-sites[0][1]) << "\n";








            //Testing for two sites additive
            vector<double> two_site_starting_params(4);
            two_site_starting_params[0] = additive_parameters[0];
            two_site_starting_params[1] = 0;
            two_site_starting_params[2] = 2 * options.site_file_morgan_positions[i] - additive_parameters[0] ;
            two_site_starting_params[3] = 0;



            sites = parameters_to_sites(two_site_starting_params, 2);


            optimizer.init_bounds(4, min(0.01, bound_size/5) );
            optimizer.min_bounds[0] = options.site_file_low_bounds[i];
            optimizer.max_bounds[0] = options.site_file_high_bounds[i];

            optimizer.min_bounds[1] = -1;

            optimizer.min_bounds[2] = options.site_file_low_bounds[i];
            optimizer.max_bounds[2] = options.site_file_high_bounds[i];

            optimizer.min_bounds[3] = -1;


            vector<double> two_site_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                two_site_bottle_necks,
                &search_sites_fast_additive,
                2
            );
            
            
            double two_site_lnl = to_be_optimized_additive(two_site_parameters);

            cout << setprecision(15) << "Two site lnl ratio: " << two_site_lnl - neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << (1-sites[0][1]) << ","<< (1-sites[0][1]/2) << ",1\t1," << (1-sites[0][1]/2)/(1-sites[0][1]) << "," << 1/(1-sites[0][1]) << "\n";
            cout << "site:\t" << sites[1][0] << "\t" << (1-sites[1][1]) << ","<< (1-sites[1][1]/2) << ",1\t1," << (1-sites[1][1]/2)/(1-sites[1][1]) << "," << 1/(1-sites[1][1]) << "\n";
            

            cerr << "Two site lnl ratio: " << two_site_lnl - neutral_lnl << "\n";
            cerr << "site:\t" << sites[0][0] << "\t" << (1-sites[0][1]) << ","<< (1-sites[0][1]/2) << ",1\t1," << (1-sites[0][1]/2)/(1-sites[0][1]) << "," << 1/(1-sites[0][1]) << "\n";
            cerr << "site:\t" << sites[1][0] << "\t" << (1-sites[1][1]) << ","<< (1-sites[1][1]/2) << ",1\t1," << (1-sites[1][1]/2)/(1-sites[1][1]) << "," << 1/(1-sites[1][1]) << "\n";
            




        





            
            cout << "Completed additive two site testing of site at " << options.site_file_positions[i] << "\n";
            cout << setprecision(15) << "Two site lnl:\t" << two_site_lnl << "\n";
            cout << setprecision(15) << "single site lnl:\t" << add_lnl << "\n";
            cout << setprecision(15) << "lnl ratio:\t" << two_site_lnl - add_lnl << "\n";

            cerr << "Completed additive two site testing of site at " << options.site_file_positions[i] << "\n";
            cerr << setprecision(15) << "Two site lnl:\t" << two_site_lnl << "\n";
            cerr << setprecision(15) << "single site lnl:\t" << add_lnl << "\n";
            cerr << setprecision(15) << "lnl ratio:\t" << two_site_lnl - add_lnl << "\n";



        }
        else if (options.site_file_options[i].compare("t") == 0){

            vector<double> two_site_starting_params(6);
            two_site_starting_params[0] = options.site_file_morgan_positions[i] - 0.01;
            two_site_starting_params[1] = 1;
            two_site_starting_params[2] = 1;
            two_site_starting_params[3] = options.site_file_morgan_positions[i] + 0.01;
            two_site_starting_params[4] = 1;
            two_site_starting_params[5] = 1;


            vector<vector<double>> sites = parameters_to_sites(two_site_starting_params);

            optimizer.init_bounds(6, min(0.01, bound_size/2) );
            optimizer.min_bounds[0] = options.site_file_low_bounds[i];
            optimizer.max_bounds[0] = options.site_file_high_bounds[i];

            optimizer.min_bounds[1] = 0;
            optimizer.min_bounds[2] = 0;

            optimizer.min_bounds[3] = options.site_file_low_bounds[i];
            optimizer.max_bounds[3] = options.site_file_high_bounds[i];

            optimizer.min_bounds[4] = 0;
            optimizer.min_bounds[5] = 0;


            vector<double> two_site_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast
            );
            
            
            double two_site_lnl = to_be_optimized(two_site_parameters);

            cerr << setprecision(15) << "Two site lnl ratio: " << two_site_lnl - neutral_lnl << "\n";
            cerr << "site:\t" << sites[0][0] << "\t" << sites[0][1] << ",1," << sites[0][2] << " or " << sites[0][1]/max(sites[0][1], sites[0][2]) << ","<< 1/max(sites[0][1], sites[0][2]) <<"," << sites[0][2]/max(sites[0][1], sites[0][2]) << "\n";
            cerr << "site:\t" << sites[1][0] << "\t" << sites[1][1] << ",1," << sites[1][2] << " or " << sites[1][1]/max(sites[1][1], sites[1][2]) << ","<< 1/max(sites[1][1], sites[1][2]) <<"," << sites[1][2]/max(sites[1][1], sites[1][2]) << "\n";

            cout << "Two site lnl ratio: " << two_site_lnl - neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << sites[0][1] << ",1," << sites[0][2] << " or " << sites[0][1]/max(sites[0][1], sites[0][2]) << ","<< 1/max(sites[0][1], sites[0][2]) <<"," << sites[0][2]/max(sites[0][1], sites[0][2]) << "\n";
            cout << "site:\t" << sites[1][0] << "\t" << sites[1][1] << ",1," << sites[1][2] << " or " << sites[1][1]/max(sites[1][1], sites[1][2]) << ","<< 1/max(sites[1][1], sites[1][2]) <<"," << sites[1][2]/max(sites[1][1], sites[1][2]) << "\n";

            

            vector<double> single_site_starting_parameters(3);
            single_site_starting_parameters[0] = options.site_file_morgan_positions[i];
            single_site_starting_parameters[1] = 1;
            single_site_starting_parameters[2] = 1;

            sites = parameters_to_sites(single_site_starting_parameters);


            optimizer.init_bounds(3, min(0.01, bound_size/2) );
            optimizer.min_bounds[0] = options.site_file_low_bounds[i];
            optimizer.max_bounds[0] = options.site_file_high_bounds[i];

            optimizer.min_bounds[1] = 0;
            optimizer.min_bounds[2] = 0;


            
            vector<double> single_site_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast
            );
            
                
            double single_site_lnl = to_be_optimized(single_site_parameters);

            cout << setprecision(15) << "\n\nSingle site lnl ratio:\t" << single_site_lnl - neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << sites[0][1] << ",1," << sites[0][2] << " or " << sites[0][1]/max(sites[0][1], sites[0][2]) << ","<< 1/max(sites[0][1], sites[0][2]) <<"," << sites[0][2]/max(sites[0][1], sites[0][2]) << "\n";

            cerr << "\n\nSingle site lnl ratio:\t" << single_site_lnl - neutral_lnl << "\n";
            cerr << "site:\t" << sites[0][0] << "\t" << sites[0][1] << ",1," << sites[0][2] << " or " << sites[0][1]/max(sites[0][1], sites[0][2]) << ","<< 1/max(sites[0][1], sites[0][2]) <<"," << sites[0][2]/max(sites[0][1], sites[0][2]) << "\n";

            
            cout << "Completed Two site testing of site at " << options.site_file_positions[i] << "\n";
            cout << setprecision(15) << "Two site lnl:\t" << two_site_lnl << "\n";
            cout << setprecision(15) << "single site lnl:\t" << single_site_lnl << "\n";

            cerr << "Completed Two site testing of site at " << options.site_file_positions[i] << "\n";
            cerr << setprecision(15) << "Two site lnl:\t" << two_site_lnl << "\n";
            cerr << setprecision(15) << "single site lnl:\t" << single_site_lnl << "\n";
        
        }else if (options.site_file_options[i].compare("S") == 0){

            //FOR TESTING PERPOSES I WILL SUPPLY THE SITES
            vector<double> loci(6);
            loci[0] = 0.05008 * 0.3025;
            loci[1] = 0.06598 * 0.3025;

            loci[2] = 0.11139 * 0.3025;
            loci[3] = 0.13544 * 0.3025;

            loci[4] = 0.27583 * 0.3025;
            loci[5] = 0.29938 * 0.3025;

            context.restricted_search_sites = loci;

            vector<double> starting_parameters(6);
            starting_parameters[0] = 0;
            starting_parameters[1] = 0;
            starting_parameters[2] = 0;
            starting_parameters[3] = 0;
            starting_parameters[4] = 0;
            starting_parameters[5] = 0;


            vector<vector<double>> sites = parameters_to_sites(starting_parameters,1);
            

            optimizer.init_bounds(6, 0.01);
            optimizer.min_bounds[0] = -1;
            optimizer.min_bounds[1] = -1;
            optimizer.min_bounds[2] = -1;
            optimizer.min_bounds[3] = -1;
            optimizer.min_bounds[4] = -1;
            optimizer.min_bounds[5] = -1;


            vector<double> found_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast_restricted_additive,
                1
            );


            cerr << "Found parameters: \n";
            cout << "Found parameters: \n";
            for(int i = 0; i < 6; i++){
                cerr << loci[i] << "\t" << found_parameters[i] << "\n";
                cout << loci[i] << "\t" << found_parameters[i] << "\n";
            }
            

            // double lnl = to_be_optimized_restricted_additive(found_parameters);


            //cerr << "Lnl: " << lnl << "\n";
            //cout << "Lnl: " << lnl << "\n";



            
            

        }

        cout << setprecision(15);
        cerr << setprecision(15);
    }

    return ans;
}

#endif
