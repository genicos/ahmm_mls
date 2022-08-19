#ifndef SIMULATE_SELECTIVE_ADMIXTURE_H
#define SIMULATE_SELECTIVE_ADMIXTURE_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <utility>
#include <cmath>
#include <string>
#include <fstream>
#include <algorithm>
#include <armadillo>

#include <pthread.h>

#define EPSILON 0.0000001


using namespace arma;
using namespace std;

double gwall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


int gamete_to_index(vector<bool> gamete){
    int ans = 0;
    for(unsigned int i = 0; i < gamete.size(); i++)
        ans += gamete[i] * (1 << i);
    
    return ans;
}

vector<bool> index_to_gamete(int index, int sites){
    vector<bool> ans(sites);
    for(int s = 0; s < sites; s++)
        ans[s] = (index >> s) % 2;
    return ans;
}





void HtoD(Row<double> *H, Row<double> *D){
    for (uint i = 0 ; i < H->n_elem; i++){
        for (uint j = 0; j < H->n_elem; j++) {
            (*D)(i*H->n_elem + j) = (*H)(i) * (*H)(j);
        }
    }
}


void normH(Row<double> *H){
    double total = 0;
    for (uint i = 0; i < H->n_elem; i++){
        total += (*H)(i);
    }
    for (uint i = 0; i < H->n_elem; i++){
        (*H)(i) = (*H)(i)/total;
    }
}


//m is proportion of ancestry type 0
vector<double> adjacent_transition_rate(vector<double> recomb_rates, vector<vector<double>> fitnesses, double m,
    double n_site1, double n_site2, int generations){
    
    
    // inserting netural sites //////////////////////////////////////////////////////////////
    int sites = fitnesses.size() + 2;                                                      //
    int n_sites[2];
    
    vector<double> neutral_fit(3);
    neutral_fit[0] = 1;
    neutral_fit[1] = 1;
    neutral_fit[2] = 1;
    
    vector<double> new_recomb_rates = recomb_rates;
    vector<vector<double>> new_fitnesses = fitnesses;

    
    //insert neutral sites into vectors
    bool inserted_first = false;
    double recomb_sum = 0;
    
    for(unsigned int i = 0; i < recomb_rates.size(); i++){
        recomb_sum += recomb_rates[i];

        if(!inserted_first){
            if(n_site1 - EPSILON <= recomb_sum){
                auto fit = new_fitnesses.begin();
                new_fitnesses.insert(fit + i, neutral_fit);
                n_sites[0] = i;

                auto rit = new_recomb_rates.begin();
                new_recomb_rates.insert(rit + i + 1, recomb_sum - n_site1);
                new_recomb_rates[i] = n_site1 - (recomb_sum - recomb_rates[i]);
                

                inserted_first = true;
            }
        }
        
        if(inserted_first){
            if(n_site2 - EPSILON <= recomb_sum){
                auto fit = new_fitnesses.begin();
                new_fitnesses.insert(fit + i + 1, neutral_fit);
                n_sites[1] = i + 1;

                auto rit = new_recomb_rates.begin();
                new_recomb_rates.insert(rit + i + 2, recomb_sum - n_site2);
                new_recomb_rates[i + 1] =  n_site2 - (recomb_sum - new_recomb_rates[i+1]);

                break;
            }
        }
    }

    fitnesses = new_fitnesses;
    recomb_rates = new_recomb_rates;                                                       //
    /////////////////////////////////////////////////////////////////////////////////////////
    
    


    int haploids = pow(2,sites);
    
    Row<double> H(haploids, fill::zeros);
    H(0) = m;
    H(haploids - 1) = 1 - m;

    mat M(haploids*haploids, haploids, fill::zeros);

    Row<double> D(haploids*haploids, fill::zeros);
    
    // populate M matrix //////////////////////////////////////////////////////////////
    //loop through all diploids                                                      //
    
    for (int i = 0; i < haploids; i++){
        for (int j = 0; j < haploids; j++){
            
            int di = i*haploids + j;
            vector<bool> chr1 = index_to_gamete(i, sites);
            vector<bool> chr2 = index_to_gamete(j, sites);

            
            double fitness = 1;

            for(int s = 0; s < sites; s++){
                if     (chr1[s] == 0 && chr2[s] == 0)
                    fitness *= fitnesses[s][0];
                else if(chr1[s] == 1 && chr2[s] == 1)
                    fitness *= fitnesses[s][2];
                else
                    fitness *= fitnesses[s][1];
            }
            
            //calculating possible recombinants
            vector<bool> recomb_chr1 = chr1;
            vector<bool> recomb_chr2 = chr2;

            // split before first site
            M(di, gamete_to_index(recomb_chr1)) += recomb_rates[0]/2 * fitness;
            M(di, gamete_to_index(recomb_chr2)) += recomb_rates[0]/2 * fitness;
            
            // split directly after site s
            for(int s = 0; s < sites; s++){
                recomb_chr1[s] = chr2[s];
                recomb_chr2[s] = chr1[s];

                M(di, gamete_to_index(recomb_chr1)) += recomb_rates[s+1]/2 * fitness;
                M(di, gamete_to_index(recomb_chr2)) += recomb_rates[s+1]/2 * fitness;

            }
            
        }
    }                                                                                //
    ///////////////////////////////////////////////////////////////////////////////////
    
    
    // Run through generations /////////////
    for(int g = 0; g < generations; g++) {//
        HtoD(&H, &D);
        H = D * M;
        normH(&H);
    }                                     //
    ////////////////////////////////////////

    
    // Calculating transition rates ///////////////////////////////
    vector<double> transitions(4);                               //
    
    for(int i = 0; i < haploids; i++){
        int anc_of_site_1 = index_to_gamete(i,sites)[n_sites[0]];
        int anc_of_site_2 = index_to_gamete(i,sites)[n_sites[1]];

        if(anc_of_site_1 == 0 && anc_of_site_2 == 0)
            transitions[0] += H(i);
        
        if(anc_of_site_1 == 0 && anc_of_site_2 == 1)
            transitions[1] += H(i);

        if(anc_of_site_1 == 1 && anc_of_site_2 == 0)
            transitions[2] += H(i);

        if(anc_of_site_1 == 1 && anc_of_site_2 == 1)
            transitions[3] += H(i);

    }                                                            //
    ///////////////////////////////////////////////////////////////
    
    
    return transitions;
}


vector<vector<double>> ancestry_trajectory(vector<double> recomb_rates, vector<vector<double>> fitnesses, double m,
    int generations){
    
    int sites = fitnesses.size();
    
    int haploids = pow(2,sites);
    
    Row<double> H(haploids, fill::zeros);
    H(0) = m;
    H(haploids - 1) = 1 - m;

    mat M(haploids*haploids, haploids, fill::zeros);

    Row<double> D(haploids*haploids, fill::zeros);

    // populate M matrix //////////////////////////////////////////////////////////////
    //loop through all diploids                                                      //
    
    for (int i = 0; i < haploids; i++){
        for (int j = 0; j < haploids; j++){

            int di = i*haploids + j;
            vector<bool> chr1 = index_to_gamete(i, sites);
            vector<bool> chr2 = index_to_gamete(j, sites);


            double fitness = 1;

            for(int s = 0; s < sites; s++){
                if     (chr1[s] == 0 && chr2[s] == 0)
                    fitness *= fitnesses[s][0];
                else if(chr1[s] == 1 && chr2[s] == 1)
                    fitness *= fitnesses[s][2];
                else
                    fitness *= fitnesses[s][1];
            }

            //calculating possible recombinants
            vector<bool> recomb_chr1 = chr1;
            vector<bool> recomb_chr2 = chr2;

            // split before first site
            M(di, gamete_to_index(recomb_chr1)) += recomb_rates[0]/2 * fitness;
            M(di, gamete_to_index(recomb_chr2)) += recomb_rates[0]/2 * fitness;

            // split directly after site s
            for(int s = 0; s < sites; s++){
                recomb_chr1[s] = chr2[s];
                recomb_chr2[s] = chr1[s];

                M(di, gamete_to_index(recomb_chr1)) += recomb_rates[s+1]/2 * fitness;
                M(di, gamete_to_index(recomb_chr2)) += recomb_rates[s+1]/2 * fitness;

            }
        }
    }                                                                                //
    ///////////////////////////////////////////////////////////////////////////////////
    
    vector<vector<double>> ancestry_trajectories(generations + 1);
    
    // Run through generations /////////////
    for(int g = 0; g < generations; g++) {//
        vector<double> this_generation_ancestry(sites);

        for(int i = 0; i < sites; i++){
            this_generation_ancestry[i] = 0;
        }

        for(int i = 0; i < haploids; i++){
            
            vector<bool> gamete = index_to_gamete(i, sites);


            for(int j = 0; j < sites; j++){
                
                
                if(!gamete[j]){
                    this_generation_ancestry[j] += H(i);
                }
            }
        }

        ancestry_trajectories[g] = this_generation_ancestry;


        HtoD(&H, &D);
        H = D * M;
        normH(&H);
    }                                     //
    ////////////////////////////////////////

    vector<double> this_generation_ancestry(sites);

    for(int i = 0; i < sites; i++){
        this_generation_ancestry[i] = 0;
    }

    for(int i = 0; i < haploids; i++){
        
        vector<bool> gamete = index_to_gamete(i, sites);


        for(int j = 0; j < sites; j++){
            
            
            if(!gamete[j]){
                this_generation_ancestry[j] += H(i);
            }
        }
    }

    ancestry_trajectories[generations] = this_generation_ancestry;

    return ancestry_trajectories;
}




















vector<double> local_ancestries;



struct intra_model_shared_info{
    long t;
    vector<double> *neutral_sites;
    vector<double> *recomb_rates_around_selected_sites;
    vector<vector<double>> *fitnesses;
    double m;
    double generations;
    vector<mat> *transition_matrices;
    int cores;
    vector<double> *selected_sites;
};

void *single_window_process(void *void_info){

    struct intra_model_shared_info *info = (struct intra_model_shared_info *)void_info;
    long t = info->t;
    
    
    for(uint i = t; i < (*info->transition_matrices).size(); i+= info->cores){

        vector<double> trans = adjacent_transition_rate(*info->recomb_rates_around_selected_sites, *info->fitnesses, info->m, (*info->neutral_sites)[i], (*info->neutral_sites)[i + 1], info->generations);

        mat transition_matrix(2,2,fill::zeros);

        transition_matrix(0,0) = trans[0]/(trans[0] + trans[1]);
        transition_matrix(0,1) = 1 - transition_matrix(0,0);
        transition_matrix(1,0) = trans[2]/(trans[2] + trans[3]);
        transition_matrix(1,1) = 1 - transition_matrix(1,0);

        (*info->transition_matrices)[i] = transition_matrix;
        local_ancestries[i] = trans[0] + trans[1];
    }

    free(void_info);
    pthread_exit(NULL);
    return NULL;
}




vector<mat> calculate_transition_rates(
  vector<double> &recomb_rates_around_neutral_sites,
  vector<double> &recomb_rates_around_selected_sites,
  vector<vector<double>> &fitnesses,
  double m,
  int generations,
  int cores
){
    vector<mat> transition_matrices(recomb_rates_around_neutral_sites.size());
    local_ancestries.clear();
    local_ancestries.resize(recomb_rates_around_neutral_sites.size());
    

    // Creating site vectors //////////////////////////////////////////////////////
    vector<double>  neutral_sites(recomb_rates_around_neutral_sites.size()  + 1);//
    vector<double> selected_sites(recomb_rates_around_selected_sites.size() - 1);

    double sum = 0;
    neutral_sites[0] = 0; //extra
    for(uint i = 0; i < recomb_rates_around_neutral_sites.size(); i++){
        sum += recomb_rates_around_neutral_sites[i];
        neutral_sites[i + 1] = sum;
    }
    sum = 0;
    for(uint i = 0; i < selected_sites.size(); i++){
        sum += recomb_rates_around_selected_sites[i];
        selected_sites[i] = sum;
    }                                                                            //
    ///////////////////////////////////////////////////////////////////////////////
    

    
    
    // Creating threads to deal with independent adjacent neutral regions //////////////////
    vector<pthread_t> threads(cores);                                                     //

    
    
    for(long t = 0; t < cores; t++){

        //Passing necessary info to thread
        struct intra_model_shared_info *this_threads_info = (struct intra_model_shared_info *)malloc(sizeof(struct intra_model_shared_info));
        this_threads_info->t = t;
        this_threads_info->neutral_sites = &neutral_sites;
        this_threads_info->recomb_rates_around_selected_sites = &recomb_rates_around_selected_sites;
        this_threads_info->fitnesses = &fitnesses;
        this_threads_info->m = m;
        this_threads_info->generations = generations;
        this_threads_info->transition_matrices = &transition_matrices;
        this_threads_info->cores = cores;

        int rc = pthread_create(&threads[t], NULL, single_window_process, (void *)this_threads_info);
        if (rc) {
            cerr << "ERROR: unable to create a thread," << rc << "\n";
            exit(-1);
        }
    }

    
    //wait for all to finish by joining them
    for (int t = 0; t < cores; t++) {
        pthread_join(threads[t], NULL);
    }                                                                                     //
    ////////////////////////////////////////////////////////////////////////////////////////
    

    
    return transition_matrices;
}






void alt_create_transition_matrix ( map<int,vector<mat> > &transition_matrix , vector<vector< map< vector<transition_information>, double > > > &transition_info, vector<double> &recombination_rate, vector<int> &positions, double &number_chromosomes, vector<mat> &transition_matrices);















double fast_transitions_radius_in_morgans = 0.05;


void *single_fast_window_process(void *void_info){

    struct intra_model_shared_info *info = (struct intra_model_shared_info *)void_info;
    long t = info->t;
    
    for(uint i = t; i < (*info->transition_matrices).size(); i+= info->cores){
        
        vector<double> pertinent_selected_sites;
        vector<vector<double>> pertinent_fitnesses;

        for(uint j = 0; j < (*info->selected_sites).size(); j++){
            double distance = (*info->neutral_sites)[i] - (*info->selected_sites)[j];
            if(distance <= fast_transitions_radius_in_morgans && distance >= -fast_transitions_radius_in_morgans){
                pertinent_selected_sites.push_back((*info->selected_sites)[j]);
                pertinent_fitnesses.push_back((*info->fitnesses)[j]);
            }
        }

        if(pertinent_selected_sites.size() >= 1){

            vector<double> pertinent_recomb_rates_around_selected_sites(pertinent_selected_sites.size() + 1);
            pertinent_recomb_rates_around_selected_sites[0] = pertinent_selected_sites[0];

            for(uint j = 1; j < pertinent_selected_sites.size(); j++){
                pertinent_recomb_rates_around_selected_sites[j] = pertinent_selected_sites[j] - pertinent_selected_sites[j - 1];
            }
            pertinent_recomb_rates_around_selected_sites[pertinent_selected_sites.size()] = 1 - pertinent_selected_sites[pertinent_selected_sites.size() - 1];


            vector<double> trans = adjacent_transition_rate(pertinent_recomb_rates_around_selected_sites, pertinent_fitnesses, info->m, (*info->neutral_sites)[i], (*info->neutral_sites)[i + 1], info->generations);

            mat transition_matrix(2,2,fill::zeros);

            transition_matrix(0,0) = trans[0]/(trans[0] + trans[1]);
            transition_matrix(0,1) = 1 - transition_matrix(0,0);
            transition_matrix(1,0) = trans[2]/(trans[2] + trans[3]);
            transition_matrix(1,1) = 1 - transition_matrix(1,0);

            (*info->transition_matrices)[i] = transition_matrix;
            local_ancestries[i] = trans[0] + trans[1];
        }else{
            mat transition_matrix(2,2,fill::zeros);

            transition_matrix(0,0) = 0.5;
            transition_matrix(0,1) = 0.5;
            transition_matrix(1,0) = 0.5;
            transition_matrix(1,1) = 0.5;

            (*info->transition_matrices)[i] = transition_matrix;
        }
    }
    
    free(void_info);
    pthread_exit(NULL);
    return NULL;
}



vector<mat> fast_transition_rates(
  vector<double> &recomb_rates_around_neutral_sites,
  vector<double> &recomb_rates_around_selected_sites,
  vector<vector<double>> &fitnesses,
  double m,
  int generations,
  int cores
) {
    vector<mat> transition_matrices(recomb_rates_around_neutral_sites.size());
    
    local_ancestries.clear();
    local_ancestries.resize(recomb_rates_around_neutral_sites.size());


    // Creating site vectors //////////////////////////////////////////////////////
    vector<double>  neutral_sites(recomb_rates_around_neutral_sites.size()  + 1);//
    vector<double> selected_sites(recomb_rates_around_selected_sites.size() - 1);

    double sum = 0;
    neutral_sites[0] = 0; //extra
    for(uint i = 0; i < recomb_rates_around_neutral_sites.size(); i++){
        sum += recomb_rates_around_neutral_sites[i];
        neutral_sites[i + 1] = sum;
    }
    sum = 0;
    for(uint i = 0; i < selected_sites.size(); i++){
        sum += recomb_rates_around_selected_sites[i];
        selected_sites[i] = sum;
    }                                                                            //
    ///////////////////////////////////////////////////////////////////////////////



    // Creating threads to deal with independent adjacent neutral regions //////////////////
    vector<pthread_t> threads(cores);

    for(long t = 0; t < cores; t++){

        //Passing necessary info to thread
        struct intra_model_shared_info *this_threads_info = (struct intra_model_shared_info *)malloc(sizeof(struct intra_model_shared_info));
        this_threads_info->t = t;
        this_threads_info->neutral_sites = &neutral_sites;
        this_threads_info->recomb_rates_around_selected_sites = &recomb_rates_around_selected_sites;
        this_threads_info->fitnesses = &fitnesses;
        this_threads_info->m = m;
        this_threads_info->generations = generations;
        this_threads_info->transition_matrices = &transition_matrices;
        this_threads_info->cores = cores;
        
        this_threads_info->selected_sites = &selected_sites;

        int rc = pthread_create(&threads[t], NULL, single_fast_window_process, (void *)this_threads_info);
        if (rc) {
            cerr << "ERROR: unable to create a thread," << rc << "\n";
            exit(-1);
        }
    }

    
    //wait for all to finish by joining them
    for (int t = 0; t < cores; t++) {
        pthread_join(threads[t], NULL);
    }                                                                                     //
    ////////////////////////////////////////////////////////////////////////////////////////

    
    
    return transition_matrices;
}



















//Actually 1 above number of pairs skipped
int pairs_skipped = 5;



//This version skips a few pairs, and interpolates in between them


void *alt_single_fast_window_process(void *void_info) {

    struct intra_model_shared_info *info = (struct intra_model_shared_info *)void_info;
    long t = info->t;
    

    //TODO each pair is responsible to spread the interpolation to those pairs around it
    
    for(uint i = t*pairs_skipped + (pairs_skipped/2); i < (*info->transition_matrices).size(); i+= info->cores*pairs_skipped){
        
        vector<double> pertinent_selected_sites;
        vector<vector<double>> pertinent_fitnesses;

        for(uint j = 0; j < (*info->selected_sites).size(); j++){
            double distance = (*info->neutral_sites)[i] - (*info->selected_sites)[j];
            if(distance <= fast_transitions_radius_in_morgans && distance >= -fast_transitions_radius_in_morgans){
                pertinent_selected_sites.push_back((*info->selected_sites)[j]);
                pertinent_fitnesses.push_back((*info->fitnesses)[j]);
            }
        }


        vector<double> pertinent_recomb_rates_around_selected_sites(pertinent_selected_sites.size() + 1);

        if(pertinent_selected_sites.size() >= 1){

            pertinent_recomb_rates_around_selected_sites[0] = pertinent_selected_sites[0];
            
            for(uint j = 1; j < pertinent_selected_sites.size(); j++){
                pertinent_recomb_rates_around_selected_sites[j] = pertinent_selected_sites[j] - pertinent_selected_sites[j - 1];
            }

            pertinent_recomb_rates_around_selected_sites[pertinent_selected_sites.size()] = 1 - pertinent_selected_sites[pertinent_selected_sites.size() - 1];
        }else{

            pertinent_recomb_rates_around_selected_sites[0] = 1;
        }

        


        vector<double> trans = adjacent_transition_rate(pertinent_recomb_rates_around_selected_sites, pertinent_fitnesses, info->m, (*info->neutral_sites)[i], (*info->neutral_sites)[i + 1], info->generations);

        mat transition_matrix(2,2,fill::zeros);

        transition_matrix(0,0) = trans[0]/(trans[0] + trans[1]);
        transition_matrix(0,1) = 1 - transition_matrix(0,0);
        transition_matrix(1,0) = trans[2]/(trans[2] + trans[3]);
        transition_matrix(1,1) = 1 - transition_matrix(1,0);


        for(int j = i - (pairs_skipped/2); j <= i + ((pairs_skipped - 1)/2) && j < local_ancestries.size(); j++){

            (*info->transition_matrices)[j] = transition_matrix;
            local_ancestries[j] = trans[0] + trans[1];
        }
        
        //If this is on the last set, then fill the rest
        if(i > local_ancestries.size() - pairs_skipped){
            for(int j = i + ((pairs_skipped - 1)/2) + 1; j < local_ancestries.size(); j++){

                (*info->transition_matrices)[j] = transition_matrix;
                local_ancestries[j] = trans[0] + trans[1];
            }
        }
        
    }
    
    free(void_info);
    pthread_exit(NULL);
    return NULL;
}


vector<mat> alternative_fast_transition_rates (
  vector<double> &recomb_rates_around_neutral_sites,
  vector<double> &recomb_rates_around_selected_sites,
  vector<vector<double>> &fitnesses,
  double m,
  int generations,
  int cores
) {
    vector<mat> transition_matrices(recomb_rates_around_neutral_sites.size());
    
    local_ancestries.clear();
    local_ancestries.resize(recomb_rates_around_neutral_sites.size());


    // Creating site vectors //////////////////////////////////////////////////////
    vector<double>  neutral_sites(recomb_rates_around_neutral_sites.size()  + 1);//
    vector<double> selected_sites(recomb_rates_around_selected_sites.size() - 1);

    double sum = 0;
    neutral_sites[0] = 0; //extra
    for(uint i = 0; i < recomb_rates_around_neutral_sites.size(); i++){
        sum += recomb_rates_around_neutral_sites[i];
        neutral_sites[i + 1] = sum;
    }
    sum = 0;
    for(uint i = 0; i < selected_sites.size(); i++){
        sum += recomb_rates_around_selected_sites[i];
        selected_sites[i] = sum;
    }                                                                            //
    ///////////////////////////////////////////////////////////////////////////////



    // Creating threads to deal with independent adjacent neutral regions //////////////////
    vector<pthread_t> threads(cores);
    
    for(long t = 0; t < cores; t++){

        //Passing necessary info to thread
        struct intra_model_shared_info *this_threads_info = (struct intra_model_shared_info *)malloc(sizeof(struct intra_model_shared_info));
        this_threads_info->t = t;
        this_threads_info->neutral_sites = &neutral_sites;
        this_threads_info->recomb_rates_around_selected_sites = &recomb_rates_around_selected_sites;
        this_threads_info->fitnesses = &fitnesses;
        this_threads_info->m = m;
        this_threads_info->generations = generations;
        this_threads_info->transition_matrices = &transition_matrices;
        this_threads_info->cores = cores;
        
        this_threads_info->selected_sites = &selected_sites;

        int rc = pthread_create(&threads[t], NULL, alt_single_fast_window_process, (void *)this_threads_info);
        if (rc) {
            cerr << "ERROR: unable to create a thread," << rc << "\n";
            exit(-1);
        }
    }

    
    //wait for all to finish by joining them
    for (int t = 0; t < cores; t++) {
        pthread_join(threads[t], NULL);
    }                                                                                     //
    ////////////////////////////////////////////////////////////////////////////////////////

    int last_present = transition_matrices.size() - 1;

    for(int i = transition_matrices.size() - 1; i >= 0; i--){
        if (transition_matrices[i].n_elem != 0){
            last_present = i;
            break;
        }
    }

    for(int i = last_present + 1; i < transition_matrices.size(); i++) {
        transition_matrices[i] = transition_matrices[last_present];
    }

    
    return transition_matrices;
}















//// create transition matrix for a given admixture model
void alt_create_transition_matrix ( map<int,vector<mat> > &transition_matrix , vector<vector< map< vector<transition_information>, double > > > &transition_info, vector<double> &recombination_rate, vector<int> &positions, double &number_chromosomes, vector<mat> &transition_matrices) {
    
    /// check if we already computed this for this sample ploidy
    if ( transition_matrix.find( number_chromosomes ) != transition_matrix.end() ) {
        return ;
    }
    
    /// else, have to create entire matrix
    /// first create data object of approporate size
    transition_matrix[number_chromosomes].resize(recombination_rate.size()) ;


    //// iterate across all positions and compute transition matrixes

    //NOTE changed p starting point from 1 to 0

    
    for ( int p = 0 ; p < recombination_rate.size(); p ++ ) {
        
        mat segment_transitions = transition_matrices[p];
        
        /// population transitions by summing across all routes
        
        transition_matrix[number_chromosomes][p] = mat(transition_info.size(),transition_info.size(), fill::zeros);

        

        for ( int i = 0 ; i < transition_info.size() ; i ++ ) {
            for ( int j = 0 ; j < transition_info[i].size() ; j ++ ) {
                
                for ( std::map<vector<transition_information>,double>::iterator t = transition_info[i][j].begin() ; t != transition_info[i][j].end() ; ++ t ) {
                    

                    double prob_t = 1 ;
                    for ( int r = 0 ; r < t->first.size() ; r ++ ) {
                        
                        prob_t *= pow( segment_transitions(t->first[r].start_state,t->first[r].end_state), t->first[r].transition_count ) ;
                        
                    }
                    
                    
                    transition_matrix[number_chromosomes][p](j,i) += prob_t * t->second ;
                }
            }
        }

    }
}



















//// create all transition rates between ancestry types for a single chromosome
mat alt_create_transition_rates ( vector<pulse> admixture_pulses, double n, vector<double> ancestry_proportion ) {
    
    /// determine ancestry proportions based on the fraction of remainder
    for ( int p = 0 ; p < admixture_pulses.size() ; p ++ ) {
        admixture_pulses[p].proportion = admixture_pulses[p].fraction_of_remainder * ancestry_proportion[admixture_pulses[p].type] ;
        ancestry_proportion[admixture_pulses[p].type] -= admixture_pulses[p].proportion ;
    }
    
    //// sort by time so oldest events are at the top
    sort( admixture_pulses.begin(), admixture_pulses.end() ) ;
    
    /// all ancestry proprionts to this point are in units of final ancestry
    /// however before the pulses that came before, there was more
    double ancestry_accounted = 0 ;
    for ( int p = admixture_pulses.size() - 1 ; p > 0 ; p -- ) {
        double accounted_here = admixture_pulses[p].proportion ;
        admixture_pulses[p].proportion = admixture_pulses[p].proportion/( 1 - ancestry_accounted ) ;
        ancestry_accounted += accounted_here ;
    }

    /// matrix to hold transition rates
    mat transition_rates( admixture_pulses.size(), admixture_pulses.size(), fill::zeros ) ;
    
    /// iterate through all states
    for ( int s1 = 0 ; s1 < admixture_pulses.size() ; s1 ++ ) {
        for ( int s2 = 0 ; s2 < admixture_pulses.size() ; s2 ++ ) {
            
            //// self transition rates are just going to be 1-all others
            if ( s2 == s1 ) continue ;
            
            /// rates calculated one way if greater than
            if ( s1 > s2 ) {
                for ( int t = s1 ; t < admixture_pulses.size() ; t ++ ) {
                    
                    /// basic recombination rates in each epoch
                    double rate ;
                    if ( t != admixture_pulses.size() - 1 ) {
                        rate = n * ( 1 - exp( (admixture_pulses[t+1].time-admixture_pulses[t].time)/n ) ) ;
                    }
                    else {
                        rate = n * ( 1 - exp( (-1*admixture_pulses[t].time)/n ) ) ;
                    }
                    
                    /// probability of no recombination in prior epochs
                    /// will skip epoch of s1 since that's in the above equation by definition
                    for ( int t2 = t - 1 ; t2 > s1 - 1 ; t2 -- ) {
                        rate *= exp((admixture_pulses[t2+1].time-admixture_pulses[t2].time)/n) ;
                    }
                    
                    /// now probability of not selecting other acnestry types
                    for ( int a = t ; a > s2 ; a -- ) {
                        rate *= ( 1 - admixture_pulses[a].proportion ) ;
                    }
                    
                    /// now select the correct ancestry type
                    if ( s2 != 0 ) {
                        rate *= admixture_pulses[s2].proportion ;
                    }
                    
                    transition_rates(admixture_pulses[s2].entry_order,admixture_pulses[s1].entry_order) += rate ;
                }
            }
            else {
                for ( int t = s2 ; t < admixture_pulses.size() ; t ++ ) {
                    
                    /// basic recombination rates in each epoch
                    double rate ;
                    if ( t != admixture_pulses.size() - 1 ) {
                        rate = n * ( 1 - exp( (admixture_pulses[t+1].time-admixture_pulses[t].time)/n ) ) ;
                    }
                    else {
                        rate = n * ( 1 - exp( (-1*admixture_pulses[t].time)/n ) ) ;
                    }
                    
                    /// probability of no recombination in prior epochs
                    /// will skip epoch of s1 since that's in the above equation by definition
                    for ( int t2 = t - 1 ; t2 > s2 - 1 ; t2 -- ) {
                        rate *= exp((admixture_pulses[t2+1].time-admixture_pulses[t2].time)/n) ;
                    }
                    
                    /// now deal with selecting lineage of the correct ancestry type
                    for ( int a = t ; a > s2 ; a -- ) {
                        rate *= ( 1 - admixture_pulses[a].proportion ) ;
                    }
                    
                    /// now select the correct ancestry type
                    rate *= admixture_pulses[s2].proportion ;
                    
                    /// and augment by this rate for this epoch
                    transition_rates(admixture_pulses[s2].entry_order,admixture_pulses[s1].entry_order) += rate ;
                }
            }
        }
    }
    
    return transition_rates.t() ;
}

#endif



//TODO
void testing_stuff(){

    //Testing that with two neutral sites,
    //the transition rates between them as calculated
    //with adjacent_transition_rate and ancestry_trajectory are the same
    
    vector<double> test_1_recomb(1);
    test_1_recomb[0] = 1;
    vector<vector<double>> test_1_fitness(0);

    int test_1_t = 10;
    double test_1_m = 0.2;


    double test_1_n1 = 0.1;
    double test_1_n2 = 0.3;
    

    vector<double> test_1_traj_recomb(3);
    test_1_traj_recomb[0] = test_1_n1;
    test_1_traj_recomb[1] = test_1_n2 - test_1_n1;
    test_1_traj_recomb[2] = 1 - test_1_n2;
    

    vector<vector<double>> test_1_traj_fit(2);
    vector<double> neut_fit(3);
    neut_fit[0] = 1;
    neut_fit[1] = 1;
    neut_fit[2] = 1;
    test_1_traj_fit[0] = neut_fit;
    test_1_traj_fit[1] = neut_fit;
    

    vector<double> test_1_adj = adjacent_transition_rate(test_1_recomb, test_1_fitness, test_1_m,
    test_1_n1, test_1_n2, test_1_t);
    

    vector<double> test_1_anc_1(2);
    test_1_anc_1[0] = test_1_adj[0] + test_1_adj[1];
    test_1_anc_1[1] = test_1_adj[0] + test_1_adj[2];

    

    vector<vector<double>> test_1_traj = ancestry_trajectory(test_1_traj_recomb, test_1_traj_fit, test_1_m,
    test_1_t);
   

    vector<double> test_1_anc_2 = test_1_traj[test_1_t - 1];







    //Testing that with a selected site,
    //the ancestry rates as calculated
    //with adjacent_transition_rate and ancestry_trajectory are the same

    int test_2_t = 265;
    double test_2_m = 0.6152;


    double test_2_n1 = 0.2354;
    double test_2_n2 = 0.2531;

    vector<double> test_2_fit(3);
    test_2_fit[0] = 0.236;
    test_2_fit[1] = 1;
    test_2_fit[2] = 1.24;


    
    vector<double> test_2_recomb(2);
    test_2_recomb[0] = test_2_n2;
    test_2_recomb[1] = 1 - test_2_n2;

    vector<vector<double>> test_2_fitness(1);
    test_2_fitness[0] = test_2_fit;

    

    vector<double> test_2_traj_recomb(3);
    test_2_traj_recomb[0] = test_2_n1;
    test_2_traj_recomb[1] = test_2_n2 - test_2_n1;
    test_2_traj_recomb[2] = 1 - test_2_n2;
   

    vector<vector<double>> test_2_traj_fit(2);
    
    test_2_traj_fit[0] = neut_fit;
    test_2_traj_fit[1] = test_2_fit;
    

    vector<double> test_2_adj = adjacent_transition_rate(test_2_recomb, test_2_fitness, test_2_m,
    test_2_n1, test_2_n2, test_2_t);
    

    vector<double> test_2_anc_1(2);
    test_2_anc_1[0] = test_2_adj[0] + test_2_adj[1];
    test_2_anc_1[1] = test_2_adj[0] + test_2_adj[2];

    

    vector<vector<double>> test_2_traj = ancestry_trajectory(test_2_traj_recomb, test_2_traj_fit, test_2_m,
    test_2_t);
    

    vector<double> test_2_anc_2 = test_2_traj[test_2_t];

}





