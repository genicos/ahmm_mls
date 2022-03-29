#ifndef MODEL_EXAMINING
#define MODEL_EXAMINING
#include <vector>
#include <math.h>

#include "optimize_selection.h"

void selection_opt::test_models(){
    double chrom_size = 0;
    for(uint i = 0; i < n_recombs.size(); i++){
        chrom_size += n_recombs[i];
    }

    context = *this;

    vector<double> empty(0);                    
    double neutral_lnl = to_be_optimized(empty);
    
    cerr << "Neutral likelihood: " << setprecision(15) << neutral_lnl << "\n";

    for(int i = 0; i < options.models.size(); i++){
        cout << setprecision(15) << to_be_optimized(options.models[i]) << "\n";
    }
}
#endif