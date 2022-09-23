#ifndef READ_SITES
#define READ_SITES
#include <string>
#include <vector>

void read_site_file(cmd_line &options, vector<double> &recomb_rates , vector<int> &positions){

    // stream in file
    ifstream in ( options.site_file.c_str() );

    vector<int> morgan_or_bp(0);

    options.site_file_positions.resize(0);
    options.site_file_options.resize(0);
    options.site_file_low_bounds.resize(0);
    options.site_file_high_bounds.resize(0);
    

    string line;
    while (getline(in, line))
    {
        istringstream iss (line);
        
        int morgan_or_bp_choice; //morgan:1   bp:0
        if (!( iss >> morgan_or_bp_choice )) { break; }
        morgan_or_bp.push_back(morgan_or_bp_choice);
        

        vector<string> searches(0);
        vector<double> loci(0);
        vector<double> low_bounds(0);
        vector<double> high_bounds(0);
        

        string option;
        if (!( iss >> option )) { break; }
        searches.push_back(option);


        double entry;

        while ( iss >> entry ) {

            loci.push_back(entry);

            if (!( iss >> option )) { 
                low_bounds.push_back(0);
                high_bounds.push_back(DBL_MAX);
                break;
            }

            //check if option is a double
            auto result = double();
            auto i = std::istringstream(option);
            i >> result;

            
            if (!i.fail() && i.eof()) {
                // if it is a double, its the bounds

                double low_bound = stod(option);
                low_bounds.push_back(low_bound);

                double high_bound;
                if (!( iss >> high_bound )) { break; }

                high_bounds.push_back(high_bound);


                if (!( iss >> option )) { break; }
                searches.push_back(option);

            }else{
                //if option is not the double, its the next option
                searches.push_back(option);

                low_bounds.push_back(-DBL_MAX);
                high_bounds.push_back(DBL_MAX);
            }
            
        }

        
        
        morgan_or_bp.push_back(morgan_or_bp_choice);

        options.site_file_options.push_back(searches);
        options.site_file_positions.push_back(loci);
        

        options.site_file_low_bounds.push_back(low_bounds);
        options.site_file_high_bounds.push_back(high_bounds);
    }



    options.site_file_morgan_positions.resize(options.site_file_positions.size());

    
    vector<double> morgan_position(recomb_rates.size());
    double sum = 0;
    for(uint i = 0; i < recomb_rates.size(); i++){
        sum += recomb_rates[i];
        morgan_position[i] = sum;
    }
    
    
    // Filling site_file_morgan_positions
    int positions_index = 0;

    for(int i = 0; i < options.site_file_positions.size() ;i++ ){
        
        if(morgan_or_bp[i] == 0) {
            
            for(int k = 0; k < options.site_file_positions[i].size(); k++) {
                
                options.site_file_morgan_positions[i] = options.site_file_positions[i];

                
                //Converting test center to morgan postion, 
                //if I cant find the bp position, exit program
                int j = 0;
                for(; j < positions.size(); j++){
                    if (positions[j] >= options.site_file_positions[i][k]) {
                        options.site_file_morgan_positions[i][k] = morgan_position[j];
                        cerr << options.site_file_morgan_positions[i][k] << "beb\n";
                        break;
                    }
                }
                if(j == positions.size()){
                    cerr << "Position not found: " << options.site_file_positions[i][k] << "\n";
                    exit(0);//TODO
                }


                //Converting low bounds to morgan postion, 
                //if I cant find the bp position, leave it as is
                j = 0;
                for(; j < positions.size(); j++){
                    if (positions[j] >= options.site_file_low_bounds[i][k]) {
                        options.site_file_low_bounds[i][k] = morgan_position[j];
                        break;
                    }
                }


                //Converting high bounds to morgan postion, 
                //if I cant find the bp position, leave it as is
                j = 0;
                for(; j < positions.size(); j++){
                    if (positions[j] >= options.site_file_high_bounds[i][k]) {
                        options.site_file_high_bounds[i][k] = morgan_position[j];
                        break;
                    }
                }
            }

        }else{
            options.site_file_morgan_positions[i] = options.site_file_positions[i];
        }

    }
    

}

#endif