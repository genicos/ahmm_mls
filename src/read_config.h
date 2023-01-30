#ifndef __READ_CONFIG_H
#define __READ_CONFIG_H
#include <string>
#include <vector>


/*
bool string_is_a_double(string par) {
    auto result = double();
    auto i = std::istringstream(par);
    i >> result;

    return !i.fail() && i.eof();
}
*/

void error_string(){
    cerr << "ERROR: malformed config file\n";
    exit(0);
}

void read_config_file(cmd_line &options, vector<double> &recomb_rates , vector<int> &positions){
    
    ifstream in ( options.site_file.c_str() );
    
    vector<int> morgan_or_bp(0);
    vector<Search> config_searches;

    string line;
    while (getline(in, line))
    {
        istringstream iss (line);
        cerr << line << "\n";

        Search this_line;

        this_line.search_string = line;

        int morgan_or_bp_choice; //morgan:1   bp:0
        if (!( iss >> morgan_or_bp_choice )) { error_string(); }
        morgan_or_bp.push_back(morgan_or_bp_choice);


        //Read in m value
        string m_string;
        if (!( iss >> m_string )) { error_string(); }
        this_line.search_m = false;
        if(m_string[0] == '(') {
            this_line.search_m = true;
            m_string = m_string.substr(1, m_string.length() - 2);
        }

        this_line.start_m = stod(m_string);


        //Read in t value
        string t_string;
        if (!( iss >> t_string )) { error_string(); }
        this_line.search_t = false;
        if(t_string[0] == '(') {
            this_line.search_t = true;
            t_string = t_string.substr(1, t_string.length() - 2);
        }

        this_line.start_t = stod(t_string);
        

        vector<string> sites_info;
        string entry;
        while(iss >> entry){
            sites_info.push_back(entry);
        }


        for(int i = 0; i < sites_info.size(); i++) {
            if (!sites_info[i].compare("l")) {
                i++;

                this_line.search_l.push_back(false);
                if(sites_info[i][0] == '(') {
                    this_line.search_l.back() = true;

                    sites_info[i] = sites_info[i].substr(1, sites_info[i].length() - 2);
                }
                this_line.start_l.push_back(stod(sites_info[i]));

                i++;

                if(i >= sites_info.size() || !sites_info[i].compare("l")){
                    this_line.min_bound_l.push_back(-1);
                    this_line.max_bound_l.push_back(-1);

                    this_line.search_h.push_back(true);
                    this_line.start_h.push_back(0.5);

                    this_line.search_s.push_back(true);
                    this_line.start_s.push_back(0);
                    this_line.s_pos.push_back(false);
                    this_line.s_neg.push_back(false);

                    i--;
                    continue;
                }

                if(string_is_a_double(sites_info[i])) {
                    this_line.min_bound_l.push_back(stod(sites_info[i]));
                    i++;
                    this_line.max_bound_l.push_back(stod(sites_info[i]));
                    i++;
                }else{
                    this_line.min_bound_l.push_back(-1);
                    this_line.max_bound_l.push_back(-1);
                }

                //cerr << "AAAA\n" << i << "\n\n";
                if(i >= sites_info.size()){
                    this_line.search_h.push_back(true);
                    this_line.start_h.push_back(0.5);

                    this_line.search_s.push_back(true);
                    this_line.start_s.push_back(0);
                    this_line.s_pos.push_back(false);
                    this_line.s_neg.push_back(false);
                    break;
                }


                if (!sites_info[i].compare("h")) {
                    i++;

                    this_line.search_h.push_back(false);
                    if(sites_info[i][0] == '(') {
                        this_line.search_h.back() = true;

                        sites_info[i] = sites_info[i].substr(1, sites_info[i].length() - 2);
                    }

                    if(sites_info[i].length() > 0)
                        this_line.start_h.push_back(stod(sites_info[i]));
                    else
                        this_line.start_h.push_back(0.5);
                    
                    i++;
                    if(i >= sites_info.size() || sites_info[i].compare("s")) {
                        this_line.search_s.push_back(true);
                        this_line.start_s.push_back(0);
                        this_line.s_pos.push_back(false);
                        this_line.s_neg.push_back(false);
                        i--;
                        continue;
                    }
                    
                    
                }else{
                    this_line.search_h.push_back(true);
                    this_line.start_h.push_back(0.5);
                }

                i++;

                this_line.search_s.push_back(false);
                this_line.s_pos.push_back(false);
                this_line.s_neg.push_back(false);

                if (sites_info[i][0] == '(') {
                    this_line.search_s.back() = true;

                    if (sites_info[i].back() != ')') {
                        if(sites_info[i].back() == '+')
                            this_line.s_pos.back() = true;
                        else
                            this_line.s_neg.back() = true;
                        
                        sites_info[i] = sites_info[i].substr(0, sites_info[i].length() - 1);
                    }
                    sites_info[i] = sites_info[i].substr(1, sites_info[i].length() - 2);
                }
                cerr << "Heres the s: " << sites_info[i] << "\n";

                if(sites_info[i].length() > 0)
                    this_line.start_s.push_back(stod(sites_info[i]));
                else
                    this_line.start_s.push_back(0);
                
            }

        }
        

        config_searches.push_back(this_line);
    }

    
    for(int i = 0; i < config_searches.size(); i++){
        
        cerr << "\n\nNEW SEARCH\n";
        cerr << "search_m: " << config_searches[i].search_m << "\n";
        cerr << "start_m:"  << config_searches[i].start_m << "\n";
        cerr << "search_t: " << config_searches[i].search_t << "\n";
        cerr << "start_t:"  << config_searches[i].start_t << "\n";
        cerr << "PRECHECK\n";
        cerr << config_searches[i].search_l.size() << "\n";
        cerr << config_searches[i].start_l.size() << "\n";
        cerr << config_searches[i].min_bound_l.size() << "\n";
        cerr << config_searches[i].max_bound_l.size() << "\n";
        cerr << config_searches[i].search_h.size() << "\n";
        cerr << config_searches[i].start_h.size() << "\n";
        cerr << config_searches[i].search_s.size() << "\n";
        cerr << config_searches[i].start_s.size() << "\n";
        cerr << config_searches[i].s_pos.size() << "\n";
        cerr << config_searches[i].s_neg.size() << "\n";

        for(int j = 0; j < config_searches[i].search_l.size(); j++){
            cerr << "\nj " << j << "\n";
            cerr << "   search_l: " << config_searches[i].search_l[j] << "\n";
            cerr << "   start_l: " << config_searches[i].start_l[j] << "\n";
            cerr << "   min_bound_l: " << config_searches[i].min_bound_l[j] << "\n";
            cerr << "   max_bound_l: " << config_searches[i].max_bound_l[j] << "\n";
            cerr << "   search_h: " << config_searches[i].search_h[j] << "\n";
            cerr << "   start_h: " << config_searches[i].start_h[j] << "\n";
            cerr << "   search_s: " << config_searches[i].search_s[j] << "\n";
            cerr << "   start_s: " << config_searches[i].start_s[j] << "\n";
            cerr << "   s_pos: " << config_searches[i].s_pos[j] << "\n";
            cerr << "   s_neg: " << config_searches[i].s_neg[j] << "\n";
        }
    }

    options.mls_searches = config_searches;

    //exit(0);
}
//s start is not present
#endif