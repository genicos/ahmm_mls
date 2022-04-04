#ifndef READ_MODEL_FILE
#define READ_MODEL_FILE
#include <string>
#include <vector>

void read_model_file(cmd_line &options){

    // stream in file
    ifstream in ( options.model_file.c_str() );

    vector<double> empty(0);

    bool first_site_of_new_model = true;
    
    while ( !in.eof() ) {
        
        string line;
        stringstream line_stream;


        getline(in, line);
        
        if (line.length() > 0){

            string site, hom0, het, hom1;
            line_stream << line;
            line_stream >> site >> hom0 >> het >> hom1;

            if ( hom0.length() == 0 || het.length() == 0 || hom1.length() == 0){
                cerr << "Model file formatted incorrectly\n";
                exit(-1);
                //TODO give more helpful error message
            }

            float het_float = stof(het);

            if(first_site_of_new_model){
                options.models.push_back(empty);
                first_site_of_new_model = false;
            }
            

            options.models[options.models.size() - 1].push_back(stof(site));
            options.models[options.models.size() - 1].push_back(stof(hom0)/het_float);
            options.models[options.models.size() - 1].push_back(stof(hom1)/het_float);
        }else{
            first_site_of_new_model = true;
        }

    }
    

}

#endif