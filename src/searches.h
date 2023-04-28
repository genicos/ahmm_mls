#ifndef __SEARCHES_H
#define __SEARCHES_H

// This object defines an MLS search
class Search{
public:
    // model file line that defines search
    string search_string;

    string search_name;

    // Decoding posterior probabilities, or printing expected local ancestry
    string posterior_options;
    
    // EPERIMENTAL, UNTESTED, may fit admixture proportion
    bool search_m;
    double start_m;

    // EPERIMENTAL, UNTESTED, may fit time since admixture
    bool search_t;
    double start_t;


    // location of selected site is fit
    vector<bool> search_l;

    // starting location of selected site
    vector<double> start_l;

    // Search bounds
    vector<double> min_bound_l;
    vector<double> max_bound_l;


    // dominance coefficient of selected site is fit
    vector<bool> search_h;

    // starting dominance coefficient value
    vector<double> start_h;
    

    // selection coefficient of selected site is fit
    vector<bool> search_s;

    // Starting selection coefficient
    vector<double> start_s;

    // Search bounds
    vector<bool> s_pos;
    vector<bool> s_neg;
};

#endif