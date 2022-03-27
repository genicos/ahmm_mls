#ifndef READ_SITES
#define READ_SITES
#include <string>
#include <vector>

void read_site_file(cmd_line &options, vector<double> &recomb_rates , vector<int> &positions){

    // stream in file
    ifstream in ( options.site_file.c_str() );

    options.site_file_positions.resize(0);
    options.site_file_morgan_positions.resize(0);

    while ( !in.eof() ) {
        int position;
        string option;
        in >> position >> option;

        options.site_file_positions.push_back(position);
        options.site_file_options.push_back(option);
    }


    //Sorting positions and corresponding options //////////////////////
    struct site_pair{                                                 //
        int A;
        string B;

        bool operator<(const site_pair x) const
            { return A < x.A;}
    };
                
    //Fill vector of structs
    vector<site_pair> site_pairs(options.site_file_positions.size());
    for(int i = 0; i < site_pairs.size(); i++){
        site_pair sp;
        sp.A = options.site_file_positions[i];
        sp.B = options.site_file_options[i];
        site_pairs[i] = sp;
    }
                
    //sort
    sort(site_pairs.begin(), site_pairs.end());
                
    //return to previous representation
    for(int i = 0; i < site_pairs.size();i++){
        options.site_file_positions[i] = site_pairs[i].A;
        options.site_file_options[i]   = site_pairs[i].B;
    }                                                                 //
    ////////////////////////////////////////////////////////////////////


    
    vector<double> morgan_position(recomb_rates.size());
    double sum = 0;
    for(uint i = 0; i < recomb_rates.size(); i++){
        sum += recomb_rates[i];
        morgan_position[i] = sum;
    }
    
    // Filling site_file_morgan_positions
    int positions_index = 0;

    for(int i = 0; i < options.site_file_positions.size() ;i++ ){
        
        int j = positions_index;
        for(; positions[j] != options.site_file_positions[i] && j < positions.size(); j++){}
        
        if(j == positions.size()){
            cerr << "Position not found: " << options.site_file_positions[i] << "\n";
        }else{
            positions_index = j;
            options.site_file_morgan_positions.push_back(morgan_position[j]);
        }
    }
    

}

#endif