#ifndef __READ_MODEL_H
#define __READ_MODEL_H
#include <string>
#include <vector>

bool string_is_a_double(string par) {
    auto result = double();
    auto i = std::istringstream(par);
    i >> result;

    return !i.fail() && i.eof();
}

void error_string(){
    cerr << "ERROR: malformed model file\n";
    exit(0);
}

void read_model_file(cmd_line &options, vector<double> &recomb_rates , vector<int> &positions) {
    
    ifstream in ( options.model_file.c_str() );
    
    vector<int> morgan_or_bp(0);
    vector<Search> config_searches;

    // read each line, each line defines an MLS search
    string line;
    while (getline(in, line))
    {
        istringstream iss (line);

        Search this_line;

        this_line.search_string = line;

        // Grab first three columns
        int morgan_or_bp_choice;
        if (!( iss >> morgan_or_bp_choice )) { error_string(); }
        morgan_or_bp.push_back(morgan_or_bp_choice);

        string search_name;
        if (!( iss >> search_name )) { error_string(); }
        this_line.search_name = search_name;

        string posterior_options; 
        if (!( iss >> posterior_options )) { error_string(); }
        this_line.posterior_options = posterior_options;


        //Read in m value
        string m_string;
        if (!( iss >> m_string )) { error_string(); }
        //Fit admixture proportion EPERIMENTAL, UNTESTED
        this_line.search_m = false;
        if(m_string[0] == '(') {
            this_line.search_m = true;
            m_string = m_string.substr(1, m_string.length() - 2);
        }
        this_line.start_m = stod(m_string);


        //Read in t value
        string t_string;
        if (!( iss >> t_string )) { error_string(); }
        //Fit time since admixture EPERIMENTAL, UNTESTED
        this_line.search_t = false;
        if(t_string[0] == '(') {
            this_line.search_t = true;
            t_string = t_string.substr(1, t_string.length() - 2);
        }
        this_line.start_t = stod(t_string);
        

        // Grab the rest of the columns, which define the selected sites in the MLS model
        vector<string> sites_info;
        string entry;
        while(iss >> entry){
            sites_info.push_back(entry);
        }

        // Parse the rest of the columns
        for(int i = 0; i < sites_info.size(); i++) {
            
            // Start of new site
            if (sites_info[i].compare("l") == 0) {
                i++;

                
                this_line.search_l.push_back(false);

                // If the location is surrounded by parenthesis, we are fitting location
                if(sites_info[i][0] == '(') {
                    this_line.search_l.back() = true;

                    sites_info[i] = sites_info[i].substr(1, sites_info[i].length() - 2);
                }
                this_line.start_l.push_back(stod(sites_info[i]));

                i++;

                if(i >= sites_info.size() || sites_info[i].compare("l") == 0){
                    // We have reached the end of defining this site

                    this_line.min_bound_l.push_back(-1);
                    this_line.max_bound_l.push_back(-1);

                    this_line.search_h.push_back(true);
                    this_line.start_h.push_back(0.5);

                    this_line.search_s.push_back(true);
                    this_line.start_s.push_back(0);
                    this_line.s_pos.push_back(false);
                    this_line.s_neg.push_back(false);

                    i--;
                    continue; // go to next site
                }

                
                if(string_is_a_double(sites_info[i])) {
                    // Bounds for location have been defined
                    this_line.min_bound_l.push_back(stod(sites_info[i]));
                    i++;
                    this_line.max_bound_l.push_back(stod(sites_info[i]));
                    i++;
                }else{
                    // No bounds defined
                    this_line.min_bound_l.push_back(-1);
                    this_line.max_bound_l.push_back(-1);
                }

                if(i >= sites_info.size()){
                    // We have reached the end of this search
                    
                    this_line.search_h.push_back(true);
                    this_line.start_h.push_back(0.5);

                    this_line.search_s.push_back(true);
                    this_line.start_s.push_back(0);
                    this_line.s_pos.push_back(false);
                    this_line.s_neg.push_back(false);

                    break;
                }
                
                
                if (sites_info[i].compare("h") == 0) {
                    // dominance coefficient is defined
                    i++;

                    this_line.search_h.push_back(false);
                    if(sites_info[i][0] == '(') { // We are fitting dominance coeffient
                        this_line.search_h.back() = true;

                        sites_info[i] = sites_info[i].substr(1, sites_info[i].length() - 2);
                    }

                    // check if starting point for h has been defined
                    if(sites_info[i].length() > 0)
                        this_line.start_h.push_back(stod(sites_info[i]));
                    else
                        this_line.start_h.push_back(0.5);
                    
                    i++;

                    if(i >= sites_info.size() || sites_info[i].compare("s")) {
                        // We have reached the end of defining this site

                        this_line.search_s.push_back(true);
                        this_line.start_s.push_back(0);
                        this_line.s_pos.push_back(false);
                        this_line.s_neg.push_back(false);
                        i--;
                        continue; // go to next site
                    }
                    
                    
                } else {
                    this_line.search_h.push_back(true);
                    this_line.start_h.push_back(0.5);
                }

                // at this point, sites_info[i] is an 's' 

                i++;
                
                this_line.search_s.push_back(false);
                this_line.s_pos.push_back(false);
                this_line.s_neg.push_back(false);

                // Check if we are searching s
                if (sites_info[i][0] == '(') {
                    this_line.search_s.back() = true;

                    if (sites_info[i].back() != ')') {
                        if(sites_info[i].back() == '+') // Restricting search to positive
                            this_line.s_pos.back() = true;
                        else                            // Restricting search to negative
                            this_line.s_neg.back() = true;
                        
                        sites_info[i] = sites_info[i].substr(0, sites_info[i].length() - 1);
                    }
                    sites_info[i] = sites_info[i].substr(1, sites_info[i].length() - 2);
                }
                
                //Check if s starting point has been defined
                if(sites_info[i].length() > 0)
                    this_line.start_s.push_back(stod(sites_info[i]));
                else
                    this_line.start_s.push_back(0);
                
            }

        }
        

        config_searches.push_back(this_line);
    }


    vector<double> morgan_position(recomb_rates.size());
    double sum = 0;
    for(uint i = 0; i < recomb_rates.size(); i++){
        sum += recomb_rates[i];
        morgan_position[i] = sum;
    }

    // Go through each search and check if sites have been defined in bp coords
    //   if so, replace bp coords with morgans
    for(int i = 0; i < config_searches.size(); i++){
        if(!morgan_or_bp[i]){

            for(int k = 0; k < config_searches[i].start_l.size(); k++) {

                //Replacing start_l with morgan position
                int j = 0;
                for(; j < positions.size(); j++){
                    if (positions[j] >= config_searches[i].start_l[k]) {
                        config_searches[i].start_l[k] = morgan_position[j];
                        break;
                    }
                }
                if(j == positions.size()){
                    cerr << "Position not found: " << config_searches[i].start_l[k] << "\n";
                    exit(0);
                }

                if (config_searches[i].min_bound_l[k] != -1) { //If min and max bound were defined

                    //Replacing min_bound_l with morgan position
                    j = 0;
                    for(; j < positions.size(); j++){
                        if (positions[j] >= config_searches[i].min_bound_l[k]) {
                            config_searches[i].min_bound_l[k] = morgan_position[j];
                            break;
                        }
                    }
                    if(j == positions.size()){
                        cerr << "Position not found: " << config_searches[i].min_bound_l[k] << "\n";
                        exit(0);
                    }

                    //Replacing max_bound_l with morgan position
                    j = 0;
                    for(; j < positions.size(); j++) {
                        if (positions[j] >= config_searches[i].max_bound_l[k]) {
                            config_searches[i].max_bound_l[k] = morgan_position[j];
                            break;
                        }
                    }
                    if(j == positions.size()) {
                        cerr << "Position not found: " << config_searches[i].max_bound_l[k] << "\n";
                        exit(0);
                    }

                }
            }
        }
    }
    
    options.mls_searches = config_searches;
}


#endif