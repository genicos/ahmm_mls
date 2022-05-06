#ifndef UNINFORMED_INFERENCE
#define UNINFORMED_INFERENCE

void selection_opt::uninformed_inference(){


    // Hyper Parameters of optimization ///////////
    double nelder_mead_reflection  = 1;
    double nelder_mead_contraction = 0.5;
    double nelder_mead_expansion   = 2;
    double nelder_mead_shrinkage   = 0.5;
    //////////////////////////////////////////////

    vector<vector<vector<double>>> bottle_necks;

    vector<vector<double>> shallow;

    vector<double> shallow_short(4);
    shallow_short[0] = 3;
    shallow_short[1] = 0.03;
    shallow_short[2] = 0.01;
    shallow_short[3] = 20;

    vector<double> shallow_medium(4);
    shallow_medium[0] = 3;
    shallow_medium[1] = 0.03;
    shallow_medium[2] = 0.03;
    shallow_medium[3] = 20;

    vector<double> shallow_tall(4);
    shallow_tall[0] = 3;
    shallow_tall[1] = 0.03;
    shallow_tall[2] = 0.05;
    shallow_tall[3] = 20;

    shallow.push_back(shallow_short);
    shallow.push_back(shallow_medium);
    shallow.push_back(shallow_tall);

    vector<vector<double>> deep;

    vector<double> deep_short(4);
    deep_short[0] = 3;
    deep_short[1] = 0.01;
    deep_short[2] = 0.005;
    deep_short[3] = 5;

    vector<double> deep_tall(4);
    deep_tall[0] = 3;
    deep_tall[1] = 0.01;
    deep_tall[2] = 0.01;
    deep_tall[3] = 5;

    deep.push_back(deep_short);
    deep.push_back(deep_tall);

    bottle_necks.push_back(shallow);
    bottle_necks.push_back(deep);



    vector<vector<vector<double>>> shallow_bottle_necks;

    vector<vector<double>> shallow1;

    //vector<double> shallow_short(4);
    shallow_short[0] = 3;
    shallow_short[1] = 0.03;
    shallow_short[2] = 0.01;
    shallow_short[3] = 20;

    //vector<double> shallow_tall(4);
    shallow_tall[0] = 3;
    shallow_tall[1] = 0.03;
    shallow_tall[2] = 0.05;
    shallow_tall[3] = 20;

    shallow1.push_back(shallow_short);
    shallow1.push_back(shallow_tall);

    vector<vector<double>> shallow2;

    //vector<double> deep_short(4);
    deep_short[0] = 3;
    deep_short[1] = 0.01;
    deep_short[2] = 0.005;
    deep_short[3] = 5;

    shallow2.push_back(deep_short);

    shallow_bottle_necks.push_back(shallow1);
    shallow_bottle_necks.push_back(shallow2);


    // Create optimizer
    nelder_mead optimizer(
        nelder_mead_reflection,
        nelder_mead_contraction,
        nelder_mead_expansion,
        nelder_mead_shrinkage
    );

    //set_context();
    context = *this;


    // Calculating lnl for neutral model
    vector<double> empty(0);
    double neutral_lnl = to_be_optimized(empty);
    

    cerr << "Neutral likelihood: " << setprecision(15) << neutral_lnl << "\n";

    neutral_transition_matrices = last_calculated_transition_matricies;

    double window_size = 0.001;



    //iterative selected site adding
    vector<vector<double>> sites;
    vector<bool> site_has_been_deep_searched; 

    double last_lnl = neutral_lnl;
    cout << "Neutral lnl\t" << setprecision(15) << neutral_lnl << "\n";
    vector<double> data_ancestry;
    vector<double> expected_ancestry;  

    data_ancestry = get_local_ancestry(last_calculated_transition_matricies);
    
    /*
    vector<double> real(3);
    real[0] = 0.2;
    real[1] = 1.00502512;
    real[2] = 0.9949748;
    double real_lnl = to_be_optimized(real);
    */



    //STANDARD: for now, only adds 5 sites
    while(sites.size() < 5){

        //data_ancestry = get_local_ancestry(last_calculated_transition_matricies);
        expected_ancestry = local_ancestries;

        vector<double> smoothed_data_ancestry(data_ancestry.size());

        double largest_deviation = 0;
        int largest_deviator = 0;
        cerr << "expected\tdata\tmorgan pos\n";

        double tot = window_size;

        for(int i = 1; i < n_recombs.size() - 1; i++){

            double total = 0;
            double count = 0;
            for(int j = i; morgan_position[i] - morgan_position[j] < 0.001 && j >= 0; j--){
                total += data_ancestry[j];
                count ++;
            }
            for(int j = i + 1; morgan_position[j] - morgan_position[i] < 0.001 && j < expected_ancestry.size(); j++){
                total += data_ancestry[j];
                count ++;
            }
            smoothed_data_ancestry[i] = total/count;


            if (abs(smoothed_data_ancestry[i] - expected_ancestry[i]) > largest_deviation){
                largest_deviation = abs(smoothed_data_ancestry[i] - expected_ancestry[i]);
                largest_deviator = i;
                //cerr << "largestdeviator" << expected_ancestry[i] << "\t" << smoothed_data_ancestry[i] << "\t" << morgan_position[i] <<"\n";
            }
            if(morgan_position[i] > tot){
                cerr << setprecision(5) << expected_ancestry[i] << "\t" << smoothed_data_ancestry[i] << "\t" << morgan_position[i] << "\n";
                tot += window_size;
            }
        }


        vector<double> new_site(3);
        new_site[0] = morgan_position[largest_deviator];
        new_site[1] = 1;
        new_site[2] = 1;
        

        cerr << "placing at :\n";
        cerr << new_site[0] << "\n";
        cout << "Placing new site at: " << new_site[0] << "\n";

        vector<double> best_parameters;
        double best_ratio = -DBL_MAX;

        vector<double> found_parameters;

        
        //generalized bottlenecking
        if(true){

            //Checking if new site is near an old site
            bool new_site_is_close_to_existing_site_that_hasnt_been_deep_searched = false;
            int site_its_close_to = 0;
            for(int s = 0; s < sites.size(); s++){
                if( abs(sites[s][0] - new_site[0]) < 0.05 && !site_has_been_deep_searched[s]){
                    new_site_is_close_to_existing_site_that_hasnt_been_deep_searched = true;
                    site_its_close_to = s;
                }
            }

            cout << "close: " << new_site_is_close_to_existing_site_that_hasnt_been_deep_searched << "\n";

            if(new_site_is_close_to_existing_site_that_hasnt_been_deep_searched){
                //pop close one from sites and place it at the end
                // then do a bottle_necks search with
                // search sites fast fix all but last

                //Remove site its close to and place it at the end
                cout << "its close\n";

                vector<double> near_site = sites[site_its_close_to];
                sites.erase(sites.begin() + site_its_close_to);
                sites.push_back(near_site);

                multi_level_optimization(
                    chrom_size,
                    optimizer,
                    sites,
                    bottle_necks,
                    &search_sites_fast_fix_all_but_last
                );

                site_has_been_deep_searched[site_has_been_deep_searched.size() - 1] = true;

            }else{
                
                cout << "its not close\n";

                sites.push_back(new_site);
                site_has_been_deep_searched.push_back(false);
                // do a bottle_necks search with
                // search sites fast fix all but last

                multi_level_optimization(
                    chrom_size,
                    optimizer,
                    sites,
                    shallow_bottle_necks,
                    &search_sites_fast_fix_all_but_last
                );
            }
            
            /*
            multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                &search_sites_fast
            );
            */


            // i need a fine tune afterward when theres more than one site

            /*
            for(int j = 0; j < bottle_necks.size(); j++){
                for(int k = 0; k < bottle_necks[j].size(); k++){
                    for(int l = 0; l < bottle_necks[j][k][0]; l++){
                        cerr << "\n SEARCH " << j << "/" << bottle_necks.size() << " " << k << "/" << bottle_necks[j].size() << " " << l << "/" << bottle_necks[j][k][0] << "\n";
                        found_parameters = search_sites_fast_fix_all_but_last(chrom_size, optimizer, sites, bottle_necks[j][k][1], bottle_necks[j][k][2], bottle_necks[j][k][3]);
                        
                        cout << "\n Result of search " << j << "/" << bottle_necks.size() << " " << k << "/" << bottle_necks[j].size() << " " << l << "/" << bottle_necks[j][k][0] << "\n";
                        cout << optimizer.max_value << "\n";
                        for(uint j = 0; j < found_parameters.size(); j++){
                            cout << found_parameters[j] << "\n";
                        }
                        cout << "\n";

                        if(optimizer.max_value > best_ratio) {
                            best_ratio = optimizer.max_value;
                            best_parameters = found_parameters;
                        }
                    }
                }
                sites = parameters_to_sites(best_parameters);
                best_ratio = -DBL_MAX;
            }
            */
        }else{

        }
        //end of else


        cerr << "\n\nBest sites so far:\n";
        for(int k = 0; k < sites.size(); k++){

            cerr << "site:\t" << sites[k][0] << "\t" << sites[k][1] << ",1," << sites[k][2] << " or " << sites[k][1]/max(sites[k][1], sites[k][2]) << ","<< 1/max(sites[k][1], sites[k][2]) <<"," << sites[k][2]/max(sites[k][1], sites[k][2]) << "\n";
        }
        cerr << "\n";

        cout << "\nnew site added:\t" << sites[sites.size() - 1][0]
        << "\n" << sites[sites.size() - 1][1]
        << ",1," << sites[sites.size() - 1][2]
        << " or\n" << sites[sites.size() - 1][1]/max(sites[sites.size() - 1][1], sites[sites.size() - 1][2])
        << ","<< 1/max(sites[sites.size() - 1][1], sites[sites.size() - 1][2])
        << "," << sites[sites.size() - 1][2]/max(sites[sites.size() - 1][1], sites[sites.size() - 1][2]) << "\n";

        best_parameters = sites_to_parameters(sites);
        double new_lnl = to_be_optimized(best_parameters);

        cout <<"lnl after new site\t" << setprecision(15) << new_lnl << "\n";
        

        cerr << "lnl: " << new_lnl << "\n";

        if(new_lnl - last_lnl < 1){
            sites.pop_back();
            break;
        }
        last_lnl = new_lnl;

        //expected_ancestry = local_ancestries;

        //exit(0);
    }

    cerr << "\n\nTERMINATED SEARCH\n\n\n";
    cout << "TERMINATED SEARCH\n";

    // I have to re-check sites by themselves to see if they really add to the lnl
        

    cerr << "Found sites:\n\n";
    for(int k = 0; k < sites.size(); k++){
        cout << "site:\t" << sites[k][0] << "\t" << sites[k][1] << ",1," << sites[k][2] << " or " << sites[k][1]/max(sites[k][1], sites[k][2]) << ","<< 1/max(sites[k][1], sites[k][2]) <<"," << sites[k][2]/max(sites[k][1], sites[k][2]) << "\n";
        cerr << "site:\t" << sites[k][0] << "\t" << sites[k][1] << ",1," << sites[k][2] << " or " << sites[k][1]/max(sites[k][1], sites[k][2]) << ","<< 1/max(sites[k][1], sites[k][2]) <<"," << sites[k][2]/max(sites[k][1], sites[k][2]) << "\n";
    }

}

#endif