#ifndef UNINFORMED_INFERENCE
#define UNINFORMED_INFERENCE

void selection_opt::uninformed_inference(){


    // Hyper Parameters of optimization ///////////
    double nelder_mead_reflection  = 1;
    double nelder_mead_contraction = 0.5;
    double nelder_mead_expansion   = 2;
    double nelder_mead_shrinkage   = 0.5;
    //////////////////////////////////////////////

    // Defining multi-level optimization parameters //
    vector<vector<vector<double>>> bottle_necks;    //

    vector<vector<double>> shallow;
    
    vector<double> shallow_short(4);
    shallow_short[0] = 3;
    shallow_short[1] = 0.01;
    shallow_short[2] = 0.005;
    shallow_short[3] = 5;

    vector<double> shallow_tall(4);
    shallow_tall[0] = 3;
    shallow_tall[1] = 0.01;
    shallow_tall[2] = 0.01;
    shallow_tall[3] = 5;

    shallow.push_back(shallow_short);
    shallow.push_back(shallow_tall);

    vector<vector<double>> deep;

    vector<double> deep_short(4);
    deep_short[0] = 3;
    deep_short[1] = 0.003;
    deep_short[2] = 0.001;
    deep_short[3] = 1;

    vector<double> deep_tall(4);
    deep_tall[0] = 3;
    deep_tall[1] = 0.003;
    deep_tall[2] = 0.005;
    deep_tall[3] = 1;

    deep.push_back(deep_short);
    deep.push_back(deep_tall);

    bottle_necks.push_back(shallow);
    bottle_necks.push_back(deep);                   //
    //////////////////////////////////////////////////


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
    double fast_neutral_lnl = to_be_optimized_only_near_sites(empty);
    
    

    cerr << "Neutral likelihood: " << setprecision(15) << neutral_lnl << "\n";

    neutral_transition_matrices = last_calculated_transition_matricies;
    context.fast_neutral_lnl = to_be_optimized_only_near_sites(empty);

    double window_size = 0.0003;



    //iterative selected site adding
    vector<vector<double>> sites;
    vector<bool> site_has_been_deep_searched; 

    double last_lnl = fast_neutral_lnl;

    cout << "Neutral lnl\t" << setprecision(15) << neutral_lnl << "\n";
    cout << "fast neutral lnl\t" << setprecision(15) << fast_neutral_lnl << "\n";

    vector<double> data_ancestry;
    vector<double> expected_ancestry;  

    data_ancestry = get_local_ancestry(last_calculated_transition_matricies);

    


































    double tot;


    
    data_ancestry = get_local_ancestry(last_calculated_transition_matricies);
    expected_ancestry = local_ancestries;

    vector<double> smoothed_data_ancestry(data_ancestry.size());
    cerr << "expected\tdata\tmorgan pos\n";

    tot = window_size;

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


        if(morgan_position[i] > tot){
            cerr << setprecision(5) << expected_ancestry[i] << "\t" << smoothed_data_ancestry[i] << "\t" << morgan_position[i] << "\n";
            cout << morgan_position[i]  << "\t" << expected_ancestry[i] << "\t" << smoothed_data_ancestry[i] <<"\n";
            tot += window_size;
        }
    }
    







    vector<double> accepted_model{
        0.10, -0.0050251256281406,
        0.11, -0.0101010101010102,
    };

    cerr << "STARING CALCULATION\n";
    double full_lnl = to_be_optimized_variations(false, false, true, false, false, false, false) (accepted_model);

    cerr << "lnl:" << full_lnl-neutral_lnl << "\n";
    cout << "lnl:" << full_lnl-neutral_lnl << "\n";
    


    data_ancestry = get_local_ancestry(last_calculated_transition_matricies);
    expected_ancestry = local_ancestries;

    cerr << "expected\tdata\tmorgan pos\n";
    cout << "expected\tdata\tmorgan pos\n";

    tot = window_size;

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


        if(morgan_position[i] > tot){
            cerr << setprecision(5) << expected_ancestry[i] << "\t" << smoothed_data_ancestry[i] << "\t" << morgan_position[i] << "\n";
            cout << morgan_position[i]  << "\t" << expected_ancestry[i] << "\t" << smoothed_data_ancestry[i] <<"\n";
            tot += window_size;
        }
    }








    /*
    vector<bool> is_loci(2);
    is_loci[0] = true;
    is_loci[1] = false;


    double best_lnl = -DBL_MAX;
    
    //So ill be given a list of sites, and ill test them iteratively
    vector<double> putative_sites(7);

    putative_sites[0] = 0.0199;
    putative_sites[1] = 0.0338;
    putative_sites[2] = 0.0409;
    putative_sites[3] = 0.0084;
    putative_sites[4] = 0.0460;
    putative_sites[5] = 0.0152;
    putative_sites[6] = 0.0513;

    /* panel01

    vector<double> putative_sites(7);

    putative_sites[0] = 0.0200;
    putative_sites[1] = 0.0338;
    putative_sites[2] = 0.0409;
    putative_sites[3] = 0.0093;
    putative_sites[4] = 0.0151;
    putative_sites[5] = 0.0460;
    putative_sites[6] = 0.0525;

    /panel 2
    vector<double> putative_sites(5);

    putative_sites[0] = 0.0182;
    putative_sites[1] = 0.0218;
    putative_sites[2] = 0.0329;
    putative_sites[3] = 0.0484;
    putative_sites[4] = 0.0607;

    

    */

    

    
    /*
    vector<double> new_site(2);
    //new_site[0] = putative_sites[0];
    //new_site[1] = 0;

    //sites.push_back(new_site);

    
    vector<double> null_model{0.0201988409220236, -0.00973395535426419};

    double cut_off = 1.7743968088179969;
    
    vector<double> alt_model_center{
        0.0182, 0,
        0.0218, 0,
    };





    cerr << "Null model\n";
    cout << "Null model\n";

    for(int i = 0; i < null_model.size()/2 ; i++){
        cerr << null_model[i*2 + 0] << ", " << null_model[i*2 + 1] << ",\n";
        cout << null_model[i*2 + 0] << ", " << null_model[i*2 + 1] << ",\n";
    }


    double null_model_lnl = to_be_optimized_variations(true, false, true, false, false) (null_model);

    cerr << "Null lnl: " << null_model_lnl - context.fast_neutral_lnl<< "\n\n";
    cout << "Null lnl: " << null_model_lnl - context.fast_neutral_lnl<< "\n\n";




    vector<vector<double>> alt_model_sites = parameters_to_sites(alt_model_center, 2);

    optimizer.init_bounds(2 * alt_model_sites.size(), 0.001);
    for ( int j = 0; j < alt_model_sites.size(); j++) {
        optimizer.min_bounds[j*2 + 0] = alt_model_sites[j][0] - 0.002;
        optimizer.max_bounds[j*2 + 0] = alt_model_sites[j][0] + 0.002;

        //optimizer.max_bounds[j*2 + 1] = 1; // for 0 and 1
        optimizer.max_bounds[j*2 + 1] = 0;
    }
    
    vector<double> alt_optimized_parameters = multi_level_optimization(
        chrom_size,
        optimizer,
        alt_model_sites,
        bottle_necks,
        to_be_optimized_variations(true, false, true, false, false),
        is_loci,
        2
    );

    cerr << "Alt  model\n";
    cout << "Alt  model\n";

    for(int i = 0; i < alt_optimized_parameters.size()/2 ; i++) {
        cerr << alt_optimized_parameters[i*2 + 0] << ", " << alt_optimized_parameters[i*2 + 1] << ",\n";
        cout << alt_optimized_parameters[i*2 + 0] << ", " << alt_optimized_parameters[i*2 + 1] << ",\n";
    }

    double alt_lnl = to_be_optimized_variations(true, false, true, false, false) (alt_optimized_parameters);

    cerr << "Alt  lnl: " << alt_lnl - context.fast_neutral_lnl<< "\n\n";
    cout << "Alt  lnl: " << alt_lnl - context.fast_neutral_lnl<< "\n\n";

    double diff = alt_lnl - null_model_lnl;

    cerr << "diff: " << diff << "\n\n";
    cout << "diff: " << diff << "\n\n";

    if(diff > cut_off){
        cerr << "ACCEPTING SITE :" << alt_model_center[alt_model_center.size() - 2] << "\n";
        cout << "ACCEPTING SITE :" << alt_model_center[alt_model_center.size() - 2] << "\n";
    }else{
        cerr << "REJECTING SITE :" << alt_model_center[alt_model_center.size() - 2] << "\n";
        cout << "REJECTING SITE :" << alt_model_center[alt_model_center.size() - 2] << "\n";
    }


    /*
    for(int i = 0; i < putative_sites.size(); i++) {

        optimizer.init_bounds(2 * sites.size(), 0.001);
        new_site[0] = putative_sites[i];
        new_site[1] = 0;

        sites.push_back(new_site);
        

        optimizer.init_bounds(2 * sites.size(), 0.001);
        for ( int j = 0; j < sites.size(); j++) {
            optimizer.min_bounds[j*2 + 0] = sites[j][0] - 0.002;
            optimizer.max_bounds[j*2 + 0] = sites[j][0] + 0.002;

            //optimizer.max_bounds[j*2 + 1] = 1; // for 0 and 1
            optimizer.max_bounds[j*2 + 1] = 0;
        }



        

        vector<vector<double>> unaltered_sites = sites;


        vector<double> singular_optimized_parameters = multi_level_optimization(
            chrom_size,
            optimizer,
            sites,
            bottle_necks,
            to_be_optimized_variations(true, false, true, false, false),
            is_loci,
            2
        );

        double this_lnl = to_be_optimized_variations(true, false, true, false, false) (singular_optimized_parameters);



        cerr << "LNL: " << this_lnl  - fast_neutral_lnl << "\n\n";
        cout << "lnl ratio: " << this_lnl  - fast_neutral_lnl  << "\n\n";


        cerr << "CURRENT SITES:\n";
        cout << "CURRENT SITES:\n";
        for ( int j = 0; j < sites.size(); j++) {
            cerr << sites[j][0] << " " << sites[j][1] << "\n";
            cout << sites[j][0] << " " << sites[j][1] << "\n";
        }
        cerr << "PREV SITES:\n";
        cout << "PREV SITES:\n";
        for ( int j = 0; j < unaltered_sites.size(); j++) {
            cerr << unaltered_sites[j][0] << " " << unaltered_sites[j][1] << "\n";
            cout << unaltered_sites[j][0] << " " << unaltered_sites[j][1] << "\n";
        }
        
        cerr << "LNL diff: " << this_lnl - best_lnl << "\n";
        cout << "LNL diff: " << this_lnl - best_lnl << "\n";
        


        if(this_lnl - best_lnl > cut_off) {
            cerr << "Adding putative site " << i << ": " << putative_sites[i] << "\n";
            cout << "Adding putative site " << i << ": " << putative_sites[i] << "\n";

            

            best_lnl = this_lnl;

            sites = parameters_to_sites(singular_optimized_parameters, 2);
            

            // at this do the iterated removal of all except last
            //

        }else{
            cerr << "dismissing putative site " << i << ": " << putative_sites[i] << "\n";
            cout << "dismissing putative site " << i << ": " << putative_sites[i] << "\n";
            unaltered_sites.pop_back();

            //sites.pop_back();
            sites = unaltered_sites;
        }


        cerr << "CURRENT SITES:\n";
        cout << "CURRENT SITES:\n";
        for ( int j = 0; j < sites.size(); j++) {
            cerr << sites[j][0] << " " << sites[j][1] << "\n";
            cout << sites[j][0] << " " << sites[j][1] << "\n";
        }
        cerr << "i: " << i << "\n";
        cout << "i: " << i << "\n";

        exit(0);

    }
    

    exit(0);
    */






















    /*

    //STANDARD: for now, only adds 5 sites
    while(sites.size() < 5){

        data_ancestry = get_local_ancestry(last_calculated_transition_matricies);
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
                cout << morgan_position[i]  << "\t" << expected_ancestry[i] << "\t" << smoothed_data_ancestry[i] <<"\n";
                tot += window_size;
            }
        }




        

        vector<double> new_site(2);
        new_site[0] = morgan_position[largest_deviator];
        new_site[1] = 0;
        

        cerr << "placing at :\n";
        cerr << new_site[0] << "\n";
        cout << "Placing new site at: " << new_site[0] << "\n";

        vector<double> best_parameters;
        double best_ratio = -DBL_MAX;

        vector<double> found_parameters;

        
    
        //Checking if new site is near an old site
        bool new_site_is_close_to_existing_site_that_hasnt_been_deep_searched = false;
        int site_its_close_to = 0;
        for(int s = 0; s < sites.size(); s++){
            if( abs(sites[s][0] - new_site[0]) < 0.002 && !site_has_been_deep_searched[s]){
                new_site_is_close_to_existing_site_that_hasnt_been_deep_searched = true;
                site_its_close_to = s;
            }
        }


        
        
        cout << "close: " << new_site_is_close_to_existing_site_that_hasnt_been_deep_searched << "\n";

        vector<vector<double>> unaltered_sites = sites;


        if(new_site_is_close_to_existing_site_that_hasnt_been_deep_searched){
            //pop close one from sites and place it at the end
            // then do a bottle_necks search with
            // search sites fast fix all but last

            //Remove site its close to and place it at the end
            cout << "its close\n";
            cerr << "its close\n";

            vector<double> near_site = sites[site_its_close_to];
            sites.erase(sites.begin() + site_its_close_to);
            sites.push_back(near_site);






            
            
            

            optimizer.init_bounds(2 * sites.size(), 0.001);
            for ( int j = 0; j < sites.size(); j++) {
                optimizer.min_bounds[j*2 + 0] = 0;
                optimizer.max_bounds[j*2 + 0] = 1;

                optimizer.max_bounds[j*2 + 1] = 1; // for 0 and 1
                //optimizer.max_bounds[j*2 + 1] = 0;
            }



            

            


            found_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                to_be_optimized_variations(true, false, true, false, false),
                is_loci,
                2
            );
            

            site_has_been_deep_searched[site_has_been_deep_searched.size() - 1] = true;

        }else{
                
            cout << "its not close\n";
            cerr << "its not close\n";

            sites.push_back(new_site);
            site_has_been_deep_searched.push_back(false);


            



            

            optimizer.init_bounds(2 * sites.size(), 0.001);
            for ( int j = 0; j < sites.size(); j++) {
                optimizer.min_bounds[j*2 + 0] = 0;
                optimizer.max_bounds[j*2 + 0] = 1;

                optimizer.max_bounds[j*2 + 1] = 1; // for 0 and 1
                //optimizer.max_bounds[j*2 + 1] = 0;
            }



            

            


            found_parameters = multi_level_optimization(
                chrom_size,
                optimizer,
                sites,
                bottle_necks,
                to_be_optimized_variations(true, false, true, false, false),
                is_loci,
                2
            );


        }

        sites = parameters_to_sites(found_parameters, 2);


        cerr << "\n\nBest sites so far:\n";
        cerr << sites.size() << "\n";
        for(int k = 0; k < sites.size(); k++){

            cerr << "site:\t" << sites[k][0] << "\t" << sites[k][1] << "\n";
        }
        cerr << "\n";

        cout << "\nnew site added:\t" << sites[sites.size() - 1][0] << " " << sites[sites.size() - 1][1] << "\n";

        best_parameters = sites_to_parameters(sites, 2);




        double new_fast_lnl = to_be_optimized_variations(true, false, true, false, false) (best_parameters);
        double new_lnl = to_be_optimized_variations(false, false, true, false, false) (best_parameters);




        cout <<"lnl after new site\t" << setprecision(15) << new_lnl << "\n";
        cout <<"fast lnl after new site\t" << setprecision(15) << new_fast_lnl << "\n";
        cout <<"fast lnl neutral ratio\t" << setprecision(15) << new_fast_lnl - fast_neutral_lnl << "\n";
        cout <<"fast lnl ratio\t" << setprecision(15) << new_fast_lnl - last_lnl << "\n";
        

        cerr << "lnl: " << new_lnl << "\n";
        cerr << "lnl: " << new_fast_lnl << "\n";

        if(new_fast_lnl - last_lnl < 3 && !new_site_is_close_to_existing_site_that_hasnt_been_deep_searched){
            sites.pop_back();
            sites = unaltered_sites;













            data_ancestry = get_local_ancestry(last_calculated_transition_matricies);
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
                    cout << morgan_position[i]  << "\t" << expected_ancestry[i] << "\t" << smoothed_data_ancestry[i] <<"\n";
                    tot += window_size;
                }
            }











            break;
        }
        last_lnl = new_fast_lnl;

        //expected_ancestry = local_ancestries;

        //exit(0);
    }

    cerr << "\n\nTERMINATED SEARCH\n\n\n";
    cout << "TERMINATED SEARCH\n";

    // I have to re-check sites by themselves to see if they really add to the lnl
        

    cerr << "Found sites:\n\n";
    for(int k = 0; k < sites.size(); k++ ) {
        cout << "site:\t" << sites[k][0] << "\t" << sites[k][1] << "\n";
        cerr << "site:\t" << sites[k][0] << "\t" << sites[k][1] << "\n";
    }
    */
    
}

#endif