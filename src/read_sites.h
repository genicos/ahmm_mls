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

    /*
    while ( !in.eof() ) {
        
        int morgan_or_bp_choice; //morgan:1   bp:0

        double position;

        string option;

        double low_end;
        double high_end;

        in >> morgan_or_bp_choice >> position >> option >> low_end >> high_end;


        morgan_or_bp.push_back(morgan_or_bp_choice);
        options.site_file_positions.push_back(position);
        options.site_file_options.push_back(option);

        options.site_file_low_bounds.push_back(low_end);
        options.site_file_high_bounds.push_back(high_end);
    }
    */
    //cerr << "BBBBBBBBBBBB\n";

    string line;
    while (getline(in, line))
    {
        cerr << line << "\n";
        istringstream iss (line);
        
        int morgan_or_bp_choice; //morgan:1   bp:0

        double position;

        string option;

        double low_end;
        double high_end;

        if (!(iss >> morgan_or_bp_choice >> position >> option >> low_end >> high_end)) { break; } // error

        
        
        morgan_or_bp.push_back(morgan_or_bp_choice);
        options.site_file_positions.push_back(position);
        options.site_file_options.push_back(option);

        options.site_file_low_bounds.push_back(low_end);
        options.site_file_high_bounds.push_back(high_end);
    }

    //cerr << "AAAAAAAAAAAA" << morgan_or_bp.size() << "\n";


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
            

            //Converting test center to morgan postion, 
            //if I cant find the bp position, exit program
            int j = 0;
            for(; j < positions.size(); j++){
                if (positions[j] >= options.site_file_positions[i]) {
                    options.site_file_morgan_positions[i] = morgan_position[j];
                    cerr << options.site_file_morgan_positions[i] << "beb\n";
                    break;
                }
            }
            if(j == positions.size()){
                cerr << "Position not found: " << options.site_file_positions[i] << "\n";
                exit(0);//TODO
            }


            //Converting low bounds to morgan postion, 
            //if I cant find the bp position, leave it as is
            j = 0;
            for(; j < positions.size(); j++){
                if (positions[j] >= options.site_file_low_bounds[i]) {
                    options.site_file_low_bounds[i] = morgan_position[j];
                    break;
                }
            }


            //Converting high bounds to morgan postion, 
            //if I cant find the bp position, leave it as is
            j = 0;
            for(; j < positions.size(); j++){
                if (positions[j] >= options.site_file_high_bounds[i]) {
                    options.site_file_high_bounds[i] = morgan_position[j];
                    break;
                }
            }

        }else{
            options.site_file_morgan_positions[i] = options.site_file_positions[i];
        }

    }

    for(int i = 0; i < options.site_file_positions.size() ;i++ ){
        if (options.site_file_low_bounds[i] < 0){
            options.site_file_low_bounds[i] = 0;
        }
    }

}

#endif