#ifndef NELDER_MEAD_GENERAL
#define NELDER_MEAD_GENERAL
#include <iostream>
#include <vector>
#include <cmath>
#include <cfloat>
#include <stdlib.h>

using namespace std;

vector<vector<double>> points;


class nelder_mead{ //This version maximizes
public:

    //Implements the nelder mead algorithm, to optimize the output of f by searching its inputs.
    // To use, ititialize an object A .
    // If you want to bound the search, call A.init_bounds(), then set A.min_bounds and A.max_bounds manually.
    // Then either set A.points manually, or call A.populate_points() . 
    // If the optimized function, f, needs additional information, then set A.supporting_info manually.
    // Call A.calculate_points() .
    // If you had set a bound and may have already broken it, call A.enforce_bounds()
    // Now you may repeatedly call A.iterate() , each time performing a single iteration of the nelder mead algorithm


    vector<vector<double>> points;   // set of points in simplex
    vector<double> point_values;     // values of simplex points in f
    double max_value = -DBL_MAX;     // of simplex
    double min_value =  DBL_MAX;     // of simplex
    double scnd_min_value = DBL_MAX; // second min value of simplex
    uint max_index = 0;              // index of point with max value
    uint min_index = 0;              // index of point with min value
    uint scnd_min_index = 0;         // index of point with second min value


    bool bounded = false;
    vector<double> max_bounds;
    vector<double> min_bounds;

    // When a point goes over a bound, it is brought back into bounds at a random place within this margin
    double random_edge_bounce = 0; 

    //number of iterations in a row in which the shrink operation was performed
    int repeated_shrinkages = 0; 

    // Nelder-Mead hyperparameters
    double reflection;
    double contraction;
    double expansion;
    double shrinkage;
    nelder_mead(double r, double c, double e, double s) {
        reflection = r;
        contraction = c;
        expansion = e;
        shrinkage = s;
    }

    nelder_mead(){
        reflection = 1;
        contraction = 0.5;
        expansion = 2;
        shrinkage = 0.5;
    }

    // Extra information needed in f
    void *supporting_info = NULL;

    void init_bounds(int dimensions, double reb);
    void enforce_bounds();

    void populate_points(int dimensions, double size, vector<double> center, vector<double> scales);
    void populate_points(int dimensions, double size, vector<double> center);

    // f is function that we which to optimize
    int calculate_points(double (*f)(vector<double>, void*));  // returns number of function calls to f
    int iterate(double (*f)(vector<double>, void*));           //


    void determine_max();
    void determine_min();
    void determine_scnd_min();


    double simplex_size();     // Largest distance between any two points in simplex
    double mean_value();       // of simplex points
    vector<double> centroid(); // of simplex
};


// Initialize boundary
// dimensions is the number of dimensions each simplex point has
// reb is random_edge_bounce, 
//   if a dimension of a point is bounded above by R, and an iteration places it above R
//   then that point will be moved to a random point in [R - random_edge_bounce, R] in that dimension
void nelder_mead::init_bounds(int dimensions, double reb){

    bounded = true;

    max_bounds.clear();
    min_bounds.clear();

    random_edge_bounce = reb;
    
    for(int i = 0; i < dimensions; i++){
        max_bounds.push_back(  DBL_MAX);
        min_bounds.push_back( -DBL_MAX);
    }
    
}

void nelder_mead::enforce_bounds(){

    for(uint i = 0; i < points.size(); i++){
        for(int d = 0; d < points[i].size(); d++) {
            if(points[i][d] >= max_bounds[d]){
                points[i][d] = max_bounds[d] - random_edge_bounce * rand() / RAND_MAX;
            }else if(points[i][d] <= min_bounds[d]){
                points[i][d] = min_bounds[d] + random_edge_bounce * rand() / RAND_MAX;
            }
        }
    }
}

//Generates a simplex
//  dimension is number of dimensions each point has
//  size is size of simplex
//  center is centroid of simplex
//  scales allows you to change how far each dimesions is from the centroid
void nelder_mead::populate_points(int dimensions, double size, vector<double> center,vector<double> scales){
    int n = dimensions;


    //Calculating simplex points///////////////////////////////////////////////
    points.clear();                                                          //

    vector<double> negative_point(n);
    vector<double> centroid_point(n);
    double negative_coord = (2 -             sqrt(4 + 4*n)) / (2 * n) * size;
    double centroid_coord = (2 - 1/(n+1.0) * sqrt(4 + 4*n)) / (2 * n) * size;
    for(int i = 0; i < n; i++){
        negative_point[i] = negative_coord;
        centroid_point[i] = centroid_coord;
    }
    points.push_back(negative_point);
    

    for(int i = 0; i < n; i++){
        vector<double> basis(n,0);
        basis[i] = size;
        points.push_back(basis);
    }                                                                        //
    ///////////////////////////////////////////////////////////////////////////


    //translating simplex around center/////////
    for(uint i = 0; i < points.size(); i++){  //
        for(int d = 0; d < n; d++){
            points[i][d] -= centroid_point[d];
        }
    }                                         //
    ////////////////////////////////////////////


    //Scaling points//////////////////////////
    for(uint i = 0; i < points.size(); i++){//
        for(int d = 0; d < n; d++){
            points[i][d] *= scales[d];
        }
    }                                       //
    //////////////////////////////////////////


    //translating points to new center////////
    for(uint i = 0; i < points.size(); i++){//
        for(int d = 0; d < n; d++){
            points[i][d] += center[d];
        }
    }                                       //
    //////////////////////////////////////////
}

void nelder_mead::populate_points(int dimensions, double size, vector<double> center){
    int n = dimensions;


    //Calculating simplex points///////////////////////////////////////////////
    points.clear();                                                          //

    vector<double> negative_point(n);
    vector<double> centroid_point(n);
    double negative_coord = (2 -             sqrt(4 + 4*n)) / (2 * n) * size;
    double centroid_coord = (2 - 1/(n+1.0) * sqrt(4 + 4*n)) / (2 * n) * size;
    for(int i = 0; i < n; i++){
        negative_point[i] = negative_coord;
        centroid_point[i] = centroid_coord;
    }
    points.push_back(negative_point);
    

    for(int i = 0; i < n; i++){
        vector<double> basis(n,0);
        basis[i] = size;
        points.push_back(basis);
    }                                                                        //
    ///////////////////////////////////////////////////////////////////////////


    //translating points to new center//////////////////////
    for(uint i = 0; i < points.size(); i++){              //
        for(int d = 0; d < n; d++){
            points[i][d] += center[d] - centroid_point[d];
        }
    }                                                     //
    ////////////////////////////////////////////////////////
}



int nelder_mead::calculate_points(double (*f)(vector<double>, void*)){
    
    point_values.resize(points.size());
    
    for(uint i = 0; i < points.size(); i++){
        point_values[i] = (*f)(points[i], supporting_info);
    }
    
    determine_max();
    determine_min();
    determine_scnd_min();
    
    repeated_shrinkages = 0;
    
    return points.size();
}


int nelder_mead::iterate(double (*f)(vector<double>, void*)){

    int function_calls = 0;

    //Calculating center of best side////////////////
    vector<double> best_centroid(points[0].size());//

    for(uint d = 0; d < points[0].size(); d++){
        for(uint i = 0; i < points.size(); i++){
            if(i != min_index){
                best_centroid[d] += points[i][d];
            }
        }
        best_centroid[d] /= points.size() - 1;
    }                                              //
    /////////////////////////////////////////////////


    //Calculating reflected point//////////////////////////////////////////////////////////////
    vector<double> reflected(points[0].size());                                              //

    for(uint d = 0; d < points[0].size(); d++){
        reflected[d] = best_centroid[d] + reflection * (best_centroid[d] - points[min_index][d]);
        
        
        if(bounded){
            if(reflected[d] >= max_bounds[d]){
                reflected[d] = max_bounds[d] - random_edge_bounce * rand() / RAND_MAX;
            }else if(reflected[d] <= min_bounds[d]){
                reflected[d] = min_bounds[d] + random_edge_bounce * rand() / RAND_MAX;
            }
        }
    }
    
    
    // calculating function value of transformed point
    double reflected_value = (*f)(reflected, supporting_info);
    function_calls++;

    if(reflected_value > scnd_min_value && reflected_value < max_value){
        points[min_index] = reflected;
        point_values[min_index] = reflected_value;

        min_value = scnd_min_value;
        min_index = scnd_min_index;

        determine_scnd_min();
        
        repeated_shrinkages = 0;
        return function_calls;
    }                                                                                        //
    ///////////////////////////////////////////////////////////////////////////////////////////


    //Calculating expanded point/////////////////////////////////////////////////////////////////
    if(reflected_value > max_value){                                                           //

        vector<double> expanded(points[0].size());
        
        for(uint d = 0; d < points[0].size(); d++){
            expanded[d] = best_centroid[d] + expansion * (reflected[d] - best_centroid[d]);

            if(bounded){
                if(expanded[d] >= max_bounds[d]){
                    expanded[d] = max_bounds[d] - random_edge_bounce * rand() / RAND_MAX;
                }else if(expanded[d] <= min_bounds[d]){
                    expanded[d] = min_bounds[d] + random_edge_bounce * rand() / RAND_MAX;
                }
            }
        }


        // calculating function value of transformed point
        double expanded_value = (*f)(expanded, supporting_info);
        function_calls++;

        repeated_shrinkages = 0;
        if(expanded_value >= reflected_value){
            points[min_index] = expanded;
            point_values[min_index] = expanded_value;

            max_value = expanded_value;
            max_index = min_index;

            min_value = scnd_min_value;
            min_index = scnd_min_index;
        
            determine_scnd_min();
            
            return function_calls;
        }else{
            points[min_index] = reflected;
            point_values[min_index] = reflected_value;

            max_value = reflected_value;
            max_index = min_index;

            min_value = scnd_min_value;
            min_index = scnd_min_index;
        
            determine_scnd_min();

            return function_calls;
        }
    }                                                                                          //
    /////////////////////////////////////////////////////////////////////////////////////////////


    //calculating contracted point///////////////////////////////////////////////////////////////
    vector<double> contracted(points[0].size());                                               //

    for(uint d = 0; d < points[0].size(); d++){
        contracted[d] = best_centroid[d] + contraction * (points[min_index][d] - best_centroid[d]);
    }


    // calculating function value of transformed point
    double contracted_value = (*f)(contracted, supporting_info);
    function_calls++;

    if(contracted_value > min_value){
        points[min_index] = contracted;
        point_values[min_index] = contracted_value;

        determine_max();
        determine_min();
        determine_scnd_min();

        repeated_shrinkages = 0;
        return function_calls;
    }                                                                                          //
    /////////////////////////////////////////////////////////////////////////////////////////////

    //Shrinking simplex toward best point/////////////////////////////////////////////////////
    for(uint i = 0; i < points.size(); i++){                                                //
        
        
        if(i != max_index){
            for(uint d = 0; d < points[0].size(); d++){
                points[i][d] = points[max_index][d] + shrinkage * (points[i][d] - points[max_index][d]);
            }


            // calculating function value of transformed point
            point_values[i] = (*f)(points[i], supporting_info);
            function_calls++;
        }
    }               
    
    repeated_shrinkages++;

    determine_max();
    determine_min();
    determine_scnd_min();                                                                   //
    //////////////////////////////////////////////////////////////////////////////////////////

    return function_calls;
}



void nelder_mead::determine_max(){
    max_value = -DBL_MAX;
    for(uint i = 0; i < points.size(); i++){
        if(point_values[i] > max_value){
            max_value = point_values[i];
            max_index = i;
        }
    }
}

void nelder_mead::determine_min(){
    min_value = DBL_MAX;
    for(uint i = 0; i < points.size(); i++){
        if(point_values[i] < min_value){
            min_value = point_values[i];
            min_index = i;
        }
    }
}


void nelder_mead::determine_scnd_min(){
    scnd_min_value = DBL_MAX;

    if(min_value == max_value){
        scnd_min_value = min_value;
        return;
    }

    for(uint i = 0; i < points.size(); i++){
        if(point_values[i] < scnd_min_value && point_values[i] > min_value){
            scnd_min_value = point_values[i];
            scnd_min_index = i;
        }
    }
}


double nelder_mead::simplex_size(){
    double ans = 0;

    for(uint i = 0; i < points.size(); i++){
        for(uint j = 0; j < points.size(); j++){
            double dist = 0;
            for(uint d = 0; d < points[0].size(); d++){
                double disp = points[i][d] - points[j][d];
                dist += disp * disp;
            }
            dist = sqrt(dist);
            ans += dist;
        }
    }

    return ans / (points.size() * (points.size() - 1));
}



double nelder_mead::mean_value(){
    double ans = 0;
    for(uint i = 0; i < points.size(); i++){
        ans += point_values[i];
    }

    return ans / points.size();
}


vector<double> nelder_mead::centroid(){
    vector<double> ans(points[0].size());

    for(uint d = 0; d < points[0].size(); d++){
        ans[d] = 0;
        for(uint i = 0; i < points.size(); i++){
            ans[d] += points[i][d];
        }
        ans[d] /= points.size();
    }

    return ans;
}

#endif