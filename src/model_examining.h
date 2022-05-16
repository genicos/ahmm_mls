#ifndef MODEL_EXAMINING
#define MODEL_EXAMINING
#include <vector>
#include <math.h>

#include "optimize_selection.h"

#include <pthread.h>

int cores;
vector<vector<double>> global_models;
vector<double> global_lnl;

void *single_thread_model_evaluation(void *thread_id){

    long t = (long)thread_id;

    //int used_cores = min(global_models.size(),cores);

    for(int i = t; i < global_models.size(); i += cores){
        
        global_lnl[i] = to_be_optimized(global_models[i]);
        
    }

    return NULL;
}



void selection_opt::test_models(){
    double chrom_size = 0;
    for(uint i = 0; i < n_recombs.size(); i++){
        chrom_size += n_recombs[i];
    }

    context = *this;

    vector<double> empty(0);
    double neutral_lnl = to_be_optimized(empty);
    
    cerr << "\nNeutral likelihood: " << setprecision(15) << neutral_lnl << "\n";

    
    cores = min((int)options.models.size(), (int)options.cores);
    global_models = options.models;
    global_lnl.resize(options.models.size());

    //set default value for lnl vector
    for(int i = 0; i < global_lnl.size(); i++){
        global_lnl[i] = 1;
    }


    
    vector<pthread_t> threads(cores);

    for(long t = 0; t < cores; t++){
        int rc = pthread_create(&threads[t], NULL, single_thread_model_evaluation, (void *)t);
        if (rc) {
            cerr << "ERROR: unable to create a thread," << rc << "\n";
            exit(-1);
        }
    }

    int lines_printed = 0;

    while(lines_printed < global_lnl.size()){
        for(int i = lines_printed; i < global_lnl.size() && global_lnl[i] != 1; i++){
            cout << setprecision(15) << global_lnl[i] << "\n";
            lines_printed++;
        }
        sleep(1);
    }
    

    //wait for all to finish by joining them
    for (int t = 0; t < cores; t++) {
        pthread_join(threads[t], NULL);
    }
    
}
#endif