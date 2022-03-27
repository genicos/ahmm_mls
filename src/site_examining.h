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
    deep_short[3] = 5;

    vector<double> deep_tall(4);
    deep_tall[0] = 5;
    deep_tall[1] = 0.01;
    deep_tall[2] = 0.01;
    deep_tall[3] = 5;

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
    double neutral_lnl = to_be_optimized(empty);
    
    cerr << "Neutral likelihood: " << setprecision(15) << neutral_lnl << "\n";

    neutral_transition_matrices = last_calculated_transition_matricies;




    for(int i = 0; i < options.site_file_morgan_positions.size(); i++){
        if (options.site_file_options[i].compare("o") == 0){

            vector<double> starting_parameters(3);
            starting_parameters[0] = options.site_file_morgan_positions[i];
            starting_parameters[1] = 1;
            starting_parameters[2] = 1;

            vector<vector<double>> sites = parameters_to_sites(starting_parameters);

            
            vector<double> best_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast
            );

                
            double lnl = to_be_optimized(best_parameters);

            cout << "\n\nUnrestricted optimization lnl:\t" << lnl << "\n";
            cout << "\n\nneutral lnl:\t" << neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << sites[0][1] << ",1," << sites[0][2] << " or " << sites[0][1]/max(sites[0][1], sites[0][2]) << ","<< 1/max(sites[0][1], sites[0][2]) <<"," << sites[0][2]/max(sites[0][1], sites[0][2]) << "\n";
            
            cerr << "\n\nUnrestricted optimization lnl:\t" << lnl << "\n";
            cerr << "\n\nneutral lnl:\t" << neutral_lnl << "\n";
            cerr << "site:\t" << sites[0][0] << "\t" << sites[0][1] << ",1," << sites[0][2] << " or " << sites[0][1]/max(sites[0][1], sites[0][2]) << ","<< 1/max(sites[0][1], sites[0][2]) <<"," << sites[0][2]/max(sites[0][1], sites[0][2]) << "\n";
        }
        else if (options.site_file_options[i].compare("d") == 0){

            //Dominance testing
            vector<double> dom_starting_parameters(2);
            dom_starting_parameters[0] = options.site_file_morgan_positions[i];
            dom_starting_parameters[1] = 1;

            vector<double> add_starting_parameters(2);
            add_starting_parameters[0] = options.site_file_morgan_positions[i];
            add_starting_parameters[1] = 0;



            // Pop0 dominant testing
            vector<vector<double>> sites = parameters_to_sites(dom_starting_parameters, 2);
                
            vector<double> pop0_dom_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast_dom0,
                2
            );

            double dom0_lnl = to_be_optimized_pop0_dominant(pop0_dom_parameters);

            cout << "\n\nPop0 dom lnl ratio: " << dom0_lnl - neutral_lnl << "\n";
            cerr << "\n\nPop0 dom lnl ratio: " << dom0_lnl - neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << "1,1," << sites[0][1] << " or " << 1/max(1.0, sites[0][1]) << "," << 1/max(1.0, sites[0][1]) <<"," << sites[0][1]/max(1.0, sites[0][1]) << "\n";
            cerr << "site:\t" << sites[0][0] << "\t" << "1,1," << sites[0][1] << " or " << 1/max(1.0, sites[0][1]) << "," << 1/max(1.0, sites[0][1]) <<"," << sites[0][1]/max(1.0, sites[0][1]) << "\n";




            // Pop1 dominant testing
            sites = parameters_to_sites(dom_starting_parameters, 2);
                
            vector<double> pop1_dom_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast_dom1,
                2
            );

            double dom1_lnl = to_be_optimized_pop1_dominant(pop1_dom_parameters);

            cout << "\n\nPop1 dom lnl ratio:\t" << dom1_lnl - neutral_lnl << "\n";
            cerr << "\n\nPop1 dom lnl ratio:\t" << dom1_lnl - neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << "1,1," << sites[0][1] << " or " << 1/max(1.0, sites[0][1]) << "," << 1/max(1.0, sites[0][1]) <<"," << sites[0][1]/max(1.0, sites[0][1]) << "\n";
            cerr << "site:\t" << sites[0][0] << "\t" << "1,1," << sites[0][1] << " or " << 1/max(1.0, sites[0][1]) << "," << 1/max(1.0, sites[0][1]) <<"," << sites[0][1]/max(1.0, sites[0][1]) << "\n";



            // Additive selection testing
            sites = parameters_to_sites(add_starting_parameters, 2);

            vector<double> additive_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast_additive,
                2
            );

            double add_lnl = to_be_optimized_additive(additive_parameters);

            cout << "\n\nAdditive lnl ratio:\t" << add_lnl - neutral_lnl << "\n";
            cerr << "\n\nAdditive lnl ratio:\t" << add_lnl - neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << (1-sites[0][1]) << ","<< (1-sites[0][1]/2) << ",1\t1," << (1-sites[0][1]/2)/(1-sites[0][1]) << "," << 1/(1-sites[0][1]) << "\n";
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

            //Testing for two sites additive
            vector<double> two_site_starting_params(4);
            two_site_starting_params[0] = options.site_file_morgan_positions[i] - 0.01;
            two_site_starting_params[1] = 0;
            two_site_starting_params[2] = options.site_file_morgan_positions[i] + 0.01;
            two_site_starting_params[3] = 0;


            vector<vector<double>> sites = parameters_to_sites(two_site_starting_params, 2);

            
            vector<double> two_site_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast_additive,
                2
            );
            
            
            double two_site_lnl = to_be_optimized_additive(two_site_parameters);

            cout << "Two site lnl ratio: " << two_site_lnl - neutral_lnl << "\n";
            cerr << "Two site lnl ratio: " << two_site_lnl - neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << (1-sites[0][1]) << ","<< (1-sites[0][1]/2) << ",1\t1," << (1-sites[0][1]/2)/(1-sites[0][1]) << "," << 1/(1-sites[0][1]) << "\n";
            cerr << "site:\t" << sites[0][0] << "\t" << (1-sites[0][1]) << ","<< (1-sites[0][1]/2) << ",1\t1," << (1-sites[0][1]/2)/(1-sites[0][1]) << "," << 1/(1-sites[0][1]) << "\n";
            cout << "site:\t" << sites[1][0] << "\t" << (1-sites[1][1]) << ","<< (1-sites[1][1]/2) << ",1\t1," << (1-sites[1][1]/2)/(1-sites[1][1]) << "," << 1/(1-sites[1][1]) << "\n";
            cerr << "site:\t" << sites[1][0] << "\t" << (1-sites[1][1]) << ","<< (1-sites[1][1]/2) << ",1\t1," << (1-sites[1][1]/2)/(1-sites[1][1]) << "," << 1/(1-sites[1][1]) << "\n";
            


            




            vector<double> add_starting_parameters(2);
            add_starting_parameters[0] = options.site_file_morgan_positions[i];
            add_starting_parameters[1] = 0;

            // single site additive selection testing
            sites = parameters_to_sites(add_starting_parameters, 2);

            vector<double> additive_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast_additive,
                2
            );

            double add_lnl = to_be_optimized_additive(additive_parameters);

            cout << "\n\nAdditive lnl ratio:\t" << add_lnl - neutral_lnl << "\n";
            cerr << "\n\nAdditive lnl ratio:\t" << add_lnl - neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << (1-sites[0][1]) << ","<< (1-sites[0][1]/2) << ",1\t1," << (1-sites[0][1]/2)/(1-sites[0][1]) << "," << 1/(1-sites[0][1]) << "\n";
            cerr << "site:\t" << sites[0][0] << "\t" << (1-sites[0][1]) << ","<< (1-sites[0][1]/2) << ",1\t1," << (1-sites[0][1]/2)/(1-sites[0][1]) << "," << 1/(1-sites[0][1]) << "\n";




            
            cout << "Completed additive two site testing of site at " << options.site_file_positions[i] << "\n";
            cout << setprecision(15) << "Two site lnl:\t" << two_site_lnl << "\n";
            cout << setprecision(15) << "single site lnl:\t" << add_lnl << "\n";

            cerr << "Completed additive two site testing of site at " << options.site_file_positions[i] << "\n";
            cerr << setprecision(15) << "Two site lnl:\t" << two_site_lnl << "\n";
            cerr << setprecision(15) << "single site lnl:\t" << add_lnl << "\n";

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

            
            vector<double> two_site_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast
            );
            
            
            double two_site_lnl = to_be_optimized(two_site_parameters);

            cerr << "Two site lnl ratio: " << two_site_lnl - neutral_lnl << "\n";
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

            
            vector<double> single_site_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast
            );

                
            double single_site_lnl = to_be_optimized(single_site_parameters);

            cout << "\n\nSingle site lnl ratio:\t" << single_site_lnl - neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << sites[0][1] << ",1," << sites[0][2] << " or " << sites[0][1]/max(sites[0][1], sites[0][2]) << ","<< 1/max(sites[0][1], sites[0][2]) <<"," << sites[0][2]/max(sites[0][1], sites[0][2]) << "\n";

            cerr << "\n\nSingle site lnl ratio:\t" << single_site_lnl - neutral_lnl << "\n";
            cerr << "site:\t" << sites[0][0] << "\t" << sites[0][1] << ",1," << sites[0][2] << " or " << sites[0][1]/max(sites[0][1], sites[0][2]) << ","<< 1/max(sites[0][1], sites[0][2]) <<"," << sites[0][2]/max(sites[0][1], sites[0][2]) << "\n";

            
            cout << "Completed Two site testing of site at " << options.site_file_positions[i] << "\n";
            cout << setprecision(15) << "Two site lnl:\t" << two_site_lnl << "\n";
            cout << setprecision(15) << "single site lnl:\t" << single_site_lnl << "\n";

            cerr << "Completed Two site testing of site at " << options.site_file_positions[i] << "\n";
            cerr << setprecision(15) << "Two site lnl:\t" << two_site_lnl << "\n";
            cerr << setprecision(15) << "single site lnl:\t" << single_site_lnl << "\n";
        }

        cout << setprecision(6);
        cerr << setprecision(6);
    }

    return ans;
}

#endif