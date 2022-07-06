#ifndef SITE_EXAMINING
#define SITE_EXAMINING
#include <vector>
#include <math.h>

#include "optimize_selection.h"



void help(double s3, double f3, mat transition_matrix1, double m, int t, double n1, double n2, double scale, bool show){
    vector<double> recomb_rates_2(2);
            recomb_rates_2[0] = s3;
            recomb_rates_2[1] = 1 - s3;

            vector<vector<double>> fitness_2(1);
            vector<double> fit_2_1(3);
            fit_2_1[0] = 1;
            fit_2_1[1] = 1 + (f3 - 1)/2;
            fit_2_1[2] = f3;
            fitness_2[0] = fit_2_1;

            vector<double> trans2 = adjacent_transition_rate(recomb_rates_2, fitness_2, m, n1, n2, t);


            mat transition_matrix2(2,2,fill::zeros);

            transition_matrix2(0,0) = trans2[0]/(trans2[0] + trans2[1]);
            transition_matrix2(0,1) = 1 - transition_matrix2(0,0);
            transition_matrix2(1,0) = trans2[2]/(trans2[2] + trans2[3]);
            transition_matrix2(1,1) = 1 - transition_matrix2(1,0);
            
            if(!show){
                cerr << setprecision(17) << f3 << "\t" << s3 << "\n";
                cerr << setprecision(17) << transition_matrix1(0,0) << "\t" << transition_matrix2(0,0) << "\t" << scale*(transition_matrix2(0,0) - transition_matrix1(0,0)) << "\n";
                cerr << setprecision(17) << transition_matrix1(1,1) << "\t" << transition_matrix2(1,1) << "\t" << scale*(transition_matrix2(1,1) - transition_matrix1(1,1)) << "\n\n";
            }else{
                double diff1 = transition_matrix2(0,0) - transition_matrix1(0,0); 
                double diff2 = transition_matrix2(1,1) - transition_matrix1(1,1);

                if(diff1 < 0 && diff2 < 0){
                    cerr << "-";
                }
                if(diff1 < 0 && diff2 > 0){
                    cerr << "/";
                }
                if(diff1 > 0 && diff2 < 0){
                    cerr << "#";
                }
                if(diff1 > 0 && diff2 > 0){
                    cerr << "+";
                }
            }
}




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
    
    cerr << "\nNeutral likelihood: " << setprecision(15) << neutral_lnl << "\n";

    neutral_transition_matrices = last_calculated_transition_matricies;


    if(false){

        //Treating two sites as one

        double s1 = 0.2;
        double s2 = 0.3;
        double f1 = 2;
        double f2 = 0.5;

        double n1 = 0.40;
        double n2 = 0.401;

        static double m = 0.20;
        int t = 300;


        double s3 = 0.265476;
        

        cerr << "A\n";
        vector<double> recomb_rates_1(3);
        recomb_rates_1[0] = s1;
        recomb_rates_1[1] = s2 - s1;
        recomb_rates_1[2] = 1 - s2;

        vector<vector<double>> fitness_1(2);
        vector<double> fit_1_1(3);
        fit_1_1[0] = 1;
        fit_1_1[1] = 1 + (f1 - 1)/2;
        fit_1_1[2] = f1;
        vector<double> fit_1_2(3);
        fit_1_2[0] = 1;
        fit_1_2[1] = 1 + (f2 - 1)/2;
        fit_1_2[2] = f2;

        fitness_1[0] = fit_1_1;
        fitness_1[1] = fit_1_2;

        vector<double> trans1 = adjacent_transition_rate(recomb_rates_1, fitness_1, m, n1, n2, t);
        

        mat transition_matrix1(2,2,fill::zeros);

        transition_matrix1(0,0) = trans1[0]/(trans1[0] + trans1[1]);
        transition_matrix1(0,1) = 1 - transition_matrix1(0,0);
        transition_matrix1(1,0) = trans1[2]/(trans1[2] + trans1[3]);
        transition_matrix1(1,1) = 1 - transition_matrix1(1,0);

        double f3 = 0.5;

        double s_start = 0.40069138142637;
        double s_step =  0.0000000000000001;
        int s_num =                     100;

        double f_start = 0.99039841958817;
        double f_step =  0.0000000000000001;
        int f_num =                     100;

        for(int i = 0; i < f_num; i++){
            if(i % 10 == 0){
                cerr << (i/10)%10;
            }else{
                cerr << " ";
            }
        }
        cerr << "\n";

        for(int i = 0; i < s_num; i++){
            s3 = s_start + i * s_step;

            for(int j = 0; j < f_num; j++){
                
                f3 = f_start + j * f_step;

                help( s3,  f3,  transition_matrix1,  m,  t,  n1,  n2, 1e8, true);
            }
            cerr << setprecision(20) << " "<< s3 <<"\n";
            
        }
        cerr << f_start << "\t" << f3 << "\n";

        s3 = 0.400691381426372;
        f3 = 0.990398419588176;

        help( s3,  f3,  transition_matrix1,  m,  t,  n1,  n2, 1, false);

        
        






        

        

       

        
        exit(0);
    }

    if(false){
        //second

        //I think it might be impossible to treat 2 sites as one

        double s4 = 0.6;
        double f4 = 1;




        double s1 = 0.2;
        double s2 = 0.3;
        double f1 = 2;
        double f2 = 0.5;

        double n1 = 0.4;
        double n2 = 0.401;

        static double m = 0.2;
        int t = 300;


        vector<double> recomb_rates_1(4);
        recomb_rates_1[0] = s1;
        recomb_rates_1[1] = s2 - s1;
        recomb_rates_1[2] = s4 - s2;
        recomb_rates_1[3] = 1 - s4;


        vector<vector<double>> fitness_1(3);

        vector<double> fit_1_1(3);
        fit_1_1[0] = 1;
        fit_1_1[1] = 1 + (f1 - 1)/2;
        fit_1_1[2] = f1;
        vector<double> fit_1_2(3);
        fit_1_2[0] = 1;
        fit_1_2[1] = 1 + (f2 - 1)/2;
        fit_1_2[2] = f2;
        vector<double> fit_1_3(3);
        fit_1_3[0] = 1;
        fit_1_3[1] = 1 + (f4 - 1)/2;
        fit_1_3[2] = f4;

        cerr << fit_1_3[2] << " hey\n";

        fitness_1[0] = fit_1_1;
        fitness_1[1] = fit_1_2;
        fitness_1[2] = fit_1_3;



        vector<double> trans1 = adjacent_transition_rate(recomb_rates_1, fitness_1, m, n1, n2, t);
        
        mat transition_matrix1(2,2,fill::zeros);

        transition_matrix1(0,0) = trans1[0]/(trans1[0] + trans1[1]);
        transition_matrix1(0,1) = 1 - transition_matrix1(0,0);
        transition_matrix1(1,0) = trans1[2]/(trans1[2] + trans1[3]);
        transition_matrix1(1,1) = 1 - transition_matrix1(1,0);


        double s3 = 0.400691381426372;
        double f3 = 0.990398419588176;


        vector<double> recomb_rates_2(3);
        recomb_rates_2[0] = s3;
        recomb_rates_2[1] = s4 - s3;
        recomb_rates_2[2] = 1 - s4;

        vector<vector<double>> fitness_2(2);
        vector<double> fit_2_1(3);
        fit_2_1[0] = 1;
        fit_2_1[1] = 1 + (f3 - 1)/2;
        fit_2_1[2] = f3;

        vector<double> fit_2_2(3);
        fit_2_2[0] = 1;
        fit_2_2[1] = 1 + (f4 - 1)/2;
        fit_2_2[2] = f4;

        fitness_2[0] = fit_2_1;
        fitness_2[1] = fit_2_2;

        vector<double> trans2 = adjacent_transition_rate(recomb_rates_2, fitness_2, m, n1, n2, t);


        mat transition_matrix2(2,2,fill::zeros);

        transition_matrix2(0,0) = trans2[0]/(trans2[0] + trans2[1]);
        transition_matrix2(0,1) = 1 - transition_matrix2(0,0);
        transition_matrix2(1,0) = trans2[2]/(trans2[2] + trans2[3]);
        transition_matrix2(1,1) = 1 - transition_matrix2(1,0);


        cerr << setprecision(17) << transition_matrix1(0,0) << "\t" << transition_matrix2(0,0) << "\t" << (transition_matrix2(0,0) - transition_matrix1(0,0)) << "\n";
        cerr << setprecision(17) << transition_matrix1(1,1) << "\t" << transition_matrix2(1,1) << "\t" << (transition_matrix2(1,1) - transition_matrix1(1,1)) << "\n\n";



        exit(0);
    }


    //TODO
    if(false){
        // stream in simulated transitions
        string file_name = "simulated_transitions";
        ifstream in ( file_name.c_str() );

        vector<mat> sim_trans(0);

        while ( !in.eof() ) {
            double t00, t01, t10, t11;

            in >> t00 >> t01 >> t10 >> t11;

            mat transition_matrix(2,2,fill::zeros);

            transition_matrix(0,0) = t00;
            transition_matrix(0,1) = t01;
            transition_matrix(1,0) = t10;
            transition_matrix(1,1) = t11;

            sim_trans.push_back(transition_matrix);
        }

        double new_lnl = compute_lnl(sim_trans);

        cerr << "neutral lnl\t" << neutral_lnl << "\n";
        cerr << "new lnl    \t" << new_lnl << "\n";
        exit(0);
    }



    //TODO
    if(true){

        vector<double> test(3);
        test[0] = 0.0501;
        test[1] = 1;
        test[2] = 1;

        vector<vector<double>> sites = parameters_to_sites(test);

            
        vector<double> best_parameters = multi_level_optimization(
            chrom_size,
            optimizer,
            sites,
            bottle_necks,
            &search_sites_fast
        );

        cerr << "DONE\nDONE\nDONE\n";
        cerr << to_be_optimized_only_near_sites(best_parameters) << "\n";

        exit(0);
    }




    for(int i = 0; i < options.site_file_morgan_positions.size(); i++){
        if (options.site_file_options[i].compare("o") == 0){

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
        else if (options.site_file_options[i].compare("d") == 0){

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

            if(options.verbose_stderr){
                cerr << "Best grid search result:\n";
                cerr << dom1_starting_parameters[0] << "\n";
                cerr << dom1_starting_parameters[1] << "\n";
            }

            sites = parameters_to_sites(dom1_starting_parameters, 2);
                
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
            cout << "site:\t" << sites[0][0] << "\t" << "1,1," << sites[0][1] << " or " << 1/max(1.0, sites[0][1]) << "," << 1/max(1.0, sites[0][1]) <<"," << sites[0][1]/max(1.0, sites[0][1]) << "\n";

            cerr << "\n\nPop1 dom lnl ratio:\t" << dom1_lnl - neutral_lnl << "\n";
            cerr << "site:\t" << sites[0][0] << "\t" << "1,1," << sites[0][1] << " or " << 1/max(1.0, sites[0][1]) << "," << 1/max(1.0, sites[0][1]) <<"," << sites[0][1]/max(1.0, sites[0][1]) << "\n";



            // Additive selection testing
            add_starting_parameters = grid_search_additive(add_starting_parameters, 0.001, 0.05, 3, 6);

            if(options.verbose_stderr){
                cerr << "Best grid search result:\n";
                cerr << add_starting_parameters[0] << "\n";
                cerr << add_starting_parameters[1] << "\n";
            }

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

            cout << setprecision(15) << "Two site lnl ratio: " << two_site_lnl - neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << (1-sites[0][1]) << ","<< (1-sites[0][1]/2) << ",1\t1," << (1-sites[0][1]/2)/(1-sites[0][1]) << "," << 1/(1-sites[0][1]) << "\n";
            cout << "site:\t" << sites[1][0] << "\t" << (1-sites[1][1]) << ","<< (1-sites[1][1]/2) << ",1\t1," << (1-sites[1][1]/2)/(1-sites[1][1]) << "," << 1/(1-sites[1][1]) << "\n";
            

            cerr << "Two site lnl ratio: " << two_site_lnl - neutral_lnl << "\n";
            cerr << "site:\t" << sites[0][0] << "\t" << (1-sites[0][1]) << ","<< (1-sites[0][1]/2) << ",1\t1," << (1-sites[0][1]/2)/(1-sites[0][1]) << "," << 1/(1-sites[0][1]) << "\n";
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

            cout << setprecision(15) << "\n\nAdditive lnl ratio:\t" << add_lnl - neutral_lnl << "\n";
            cout << "site:\t" << sites[0][0] << "\t" << (1-sites[0][1]) << ","<< (1-sites[0][1]/2) << ",1\t1," << (1-sites[0][1]/2)/(1-sites[0][1]) << "," << 1/(1-sites[0][1]) << "\n";

            cerr << "\n\nAdditive lnl ratio:\t" << add_lnl - neutral_lnl << "\n";
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
        }

        cout << setprecision(6);
        cerr << setprecision(6);
    }

    return ans;
}

#endif
