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

#include <sys/time.h>

using namespace arma;
using namespace std;


double gwall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}


// Covert gamete to index in haploid vector
int gamete_to_index(vector<bool> gamete){
    int ans = 0;
    for(unsigned int i = 0; i < gamete.size(); i++)
        ans += gamete[i] * (1 << i);
    
    return ans;
}

// Covert index in haploid vector to gamete
vector<bool> index_to_gamete(int index, int sites){
    vector<bool> ans(sites);
    for(int s = 0; s < sites; s++)
        ans[s] = (index >> s) % 2;
    return ans;
}


// Fertilization, produces distribution of diploypes D as the zygotes, 
//  from distribution of haplotypes H acting as gametes. Assumes mendelian segregation
void HtoD(Row<double> *H, Row<double> *D){
    for (uint i = 0 ; i < H->n_elem; i++){
        for (uint j = 0; j < H->n_elem; j++) {
            (*D)(i*H->n_elem + j) = (*H)(i) * (*H)(j);
        }
    }
}


// Normalize haploid distribution
void normH(Row<double> *H){
    double total = 0;
    for (uint i = 0; i < H->n_elem; i++){
        total += (*H)(i);
    }
    for (uint i = 0; i < H->n_elem; i++){
        (*H)(i) = (*H)(i)/total;
    }
}


//Calculates frequency of combinations of ancestry types of adjacent sampled sites
// Returns a vector of size 4, whose entries correspond to the relative frequencies
//  of ancestry types (0,0) (0,1) (1,0) (1,1). 
//
// recomb_rates is recomb rates between selected sites
// fitnesses are relative fitnesses of genotypes of selected sites
// m is proportion of final ancestry that was introgressed from ancestral population 0
// n_site1 and n_site2 are morgan positions of adjacent sampled sites
// generations is time since admixture
vector<double> adjacent_transition_rate(vector<double> recomb_rates, vector<vector<double>> fitnesses, double m,
    double n_site1, double n_site2, int generations) {
    
    
    // inserting netural sites //////////////////////////////////////////////////////////////
    int sites = fitnesses.size() + 2;                                                      //
    int n_sites[2];
    
    vector<double> neutral_fit(3); // Fitness of neutral sites
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

        if(!inserted_first){ // We have yet to insert first neutral site
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
        
        if(inserted_first){ // We have inserted first neutral site, now inserting second
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
    
    // haploid distribution at generation 0, 
    //  only haploid are those that consist only of a single ancestry type
    Row<double> H(haploids, fill::zeros);
    H(0) = m;
    H(haploids - 1) = 1 - m;

    mat M(haploids*haploids, haploids, fill::zeros);

    Row<double> D(haploids*haploids, fill::zeros);
    
    // populate M matrix //////////////////////////////////////////////////////////////
    //loop through all diplotype                                                     //
    
    for (int i = 0; i < haploids; i++){
        for (int j = 0; j < haploids; j++){
            
            int di = i*haploids + j; // index of diplotype
            vector<bool> chr1 = index_to_gamete(i, sites); // First chromosome of diplotype
            vector<bool> chr2 = index_to_gamete(j, sites); // Second chromosome of diplotype

            
            double fitness = 1;

            // Calcuating fitness of diplotype, given relative fitness of all selected sites
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
        HtoD(&H, &D);   // fertilization
        H = D * M;      // life of diploids and meiosis
        normH(&H);
    }                                     //
    ////////////////////////////////////////

    
    // Calculating sampled site haploid frequencies ///////////////
    vector<double> ancestry_pair_frequencies(4);                 //
    
    for(int i = 0; i < haploids; i++){
        int anc_of_site_1 = index_to_gamete(i,sites)[n_sites[0]];
        int anc_of_site_2 = index_to_gamete(i,sites)[n_sites[1]];

        if(anc_of_site_1 == 0 && anc_of_site_2 == 0)
            ancestry_pair_frequencies[0] += H(i);
        
        if(anc_of_site_1 == 0 && anc_of_site_2 == 1)
            ancestry_pair_frequencies[1] += H(i);

        if(anc_of_site_1 == 1 && anc_of_site_2 == 0)
            ancestry_pair_frequencies[2] += H(i);

        if(anc_of_site_1 == 1 && anc_of_site_2 == 1)
            ancestry_pair_frequencies[3] += H(i);

    }                                                            //
    ///////////////////////////////////////////////////////////////
    
    
    return ancestry_pair_frequencies;
}















// Expected local ancestry of last calculated MLS model
vector<double> local_ancestries;


// Information needed by each thread
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
    double fast_transitions_radius_in_morgans;
    int pairs_skipped;
    double chrom_size;
};

// Each thread calls this function
void *single_window_process(void *void_info) {

    struct intra_model_shared_info *info = (struct intra_model_shared_info *)void_info;
    long t = info->t;
    
    for(uint i = t; i < (*info->transition_matrices).size(); i+= info->cores) {

        vector<double> trans = adjacent_transition_rate(*info->recomb_rates_around_selected_sites, *info->fitnesses, info->m, (*info->neutral_sites)[i], (*info->neutral_sites)[i + 1], info->generations);

        mat transition_matrix(2,2,fill::zeros);

        // Convert ancestry frequencies into transition matrix
        transition_matrix(0,0) = trans[0]/(trans[0] + trans[1]);
        transition_matrix(0,1) = 1 - transition_matrix(0,0);
        transition_matrix(1,0) = trans[2]/(trans[2] + trans[3]);
        transition_matrix(1,1) = 1 - transition_matrix(1,0);

        (*info->transition_matrices)[i] = transition_matrix;
        local_ancestries[i] = trans[0] + trans[1];

        //if verbose
        if(t == 0) {
            cerr << "thread " << t << " is at " << i << "/" << (*info->transition_matrices).size() << "\n";
        }
        
    }

    free(void_info);
    pthread_exit(NULL);
    return NULL;
}



// Calculates transition rates between all sites
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
    for(uint i = 0; i < recomb_rates_around_neutral_sites.size(); i++) {
        sum += recomb_rates_around_neutral_sites[i];
        neutral_sites[i + 1] = sum;
    }
    sum = 0;
    for(uint i = 0; i < selected_sites.size(); i++) {
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







//This version skips a few pairs, and interpolates in between them
void *alt_single_fast_window_process(void *void_info) {

    struct intra_model_shared_info *info = (struct intra_model_shared_info *)void_info;
    long t = info->t;

    int pairs_skipped = info->pairs_skipped + 1;
    
    
    for(uint i = t*pairs_skipped + (pairs_skipped/2); i < (*info->transition_matrices).size(); i+= info->cores*pairs_skipped ) {
        
        vector<double> pertinent_selected_sites;
        vector<vector<double>> pertinent_fitnesses;

        // Check which sites we are going to account for
        for(uint j = 0; j < (*info->selected_sites).size(); j++){
            double distance = (*info->neutral_sites)[i] - (*info->selected_sites)[j];
            if(distance <= info->fast_transitions_radius_in_morgans && distance >= -info->fast_transitions_radius_in_morgans){
                pertinent_selected_sites.push_back((*info->selected_sites)[j]);
                pertinent_fitnesses.push_back((*info->fitnesses)[j]);
            }
        }


        // Recomb rates around sites we are accounting for
        vector<double> pertinent_recomb_rates_around_selected_sites(pertinent_selected_sites.size() + 1);

        if(pertinent_selected_sites.size() >= 1){

            pertinent_recomb_rates_around_selected_sites[0] = pertinent_selected_sites[0];
            
            for(uint j = 1; j < pertinent_selected_sites.size(); j++){
                pertinent_recomb_rates_around_selected_sites[j] = pertinent_selected_sites[j] - pertinent_selected_sites[j - 1];
            }

            pertinent_recomb_rates_around_selected_sites[pertinent_selected_sites.size()] = info->chrom_size - pertinent_selected_sites[pertinent_selected_sites.size() - 1];
        }else{

            pertinent_recomb_rates_around_selected_sites[0] = info->chrom_size;
        }


        
        vector<double> trans = adjacent_transition_rate(pertinent_recomb_rates_around_selected_sites, pertinent_fitnesses, info->m, (*info->neutral_sites)[i], (*info->neutral_sites)[i + 1], info->generations);

        mat transition_matrix(2,2,fill::zeros);

        // Convert ancestry frequencies into transition matrix
        transition_matrix(0,0) = trans[0]/(trans[0] + trans[1]);
        transition_matrix(0,1) = 1 - transition_matrix(0,0);
        transition_matrix(1,0) = trans[2]/(trans[2] + trans[3]);
        transition_matrix(1,1) = 1 - transition_matrix(1,0);
        

        for(int j = i - (pairs_skipped/2); j <= i + ((pairs_skipped - 1)/2) && j < local_ancestries.size(); j++) {

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

// Calculates transition rates between adjacent sampled sites in a faster way
vector<mat> alternative_fast_transition_rates (
  vector<double> &recomb_rates_around_neutral_sites,
  vector<double> &recomb_rates_around_selected_sites,
  vector<vector<double>> &fitnesses,
  double m,
  int generations,
  int cores,
  double chrom_size,
  double fast_transitions_radius_in_morgans,
  int pairs_skipped
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
        this_threads_info->fast_transitions_radius_in_morgans = fast_transitions_radius_in_morgans;
        this_threads_info->pairs_skipped = pairs_skipped;
        this_threads_info->chrom_size = chrom_size;
        
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

    // Find the last transition matrix we calculated
    for(int i = transition_matrices.size() - 1; i >= 0; i--){
        if (transition_matrices[i].n_elem != 0){
            last_present = i;
            break;
        }
    }

    // Fill up the rest of the empty transition matricies with the last calculated one
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


#endif
