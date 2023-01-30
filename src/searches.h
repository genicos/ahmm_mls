#ifndef __SEARCHES_H
#define __SEARCHES_H

class Search{
public:
    string search_string;
    
    bool search_m;
    double start_m;

    bool search_t;
    double start_t;

    vector<bool> search_l;
    vector<double> start_l;
    vector<double> min_bound_l;
    vector<double> max_bound_l;

    vector<bool> search_h;
    vector<double> start_h;
    
    vector<bool> search_s;
    vector<double> start_s;
    vector<bool> s_pos;
    vector<bool> s_neg;
};

#endif