#ifndef GRID_SEARCH
#define GRID_SEARCH
#include <vector>
#include <math.h>

vector<double> selection_opt::grid_search(vector<double> center_point, double width, double height, int x_ticks, int y_ticks){
    
    vector<double> ans(3);
    double best_lnl = -DBL_MAX;

    for(int i = 0; i < x_ticks; i++){
        double site = (2.0 * i / (x_ticks-1) - 1) * abs(2.0 * i / (x_ticks-1) - 1) * width + center_point[0];
        for(int j = 0; j < y_ticks; j++){
            double f1 = (2.0 * j / (y_ticks-1) - 1) * abs(2.0 * j / (y_ticks-1) - 1) * height + center_point[1];
            for(int k = 0; k < y_ticks; k++){
                double f2 = (2.0 * k / (y_ticks-1) - 1) * abs(2.0 * k / (y_ticks-1) - 1) * height + center_point[2];
                vector<double> point(3);
                point[0] = site;
                point[1] = f1;
                point[2] = f2;
                
                double this_lnl = to_be_optimized_only_near_sites(point);
                if(this_lnl > best_lnl){
                    best_lnl = this_lnl;
                    ans = point;
                }
            }
        }
    }

    return ans;
}

vector<double> selection_opt::grid_search_additive(vector<double> center_point, double width, double height, int x_ticks, int y_ticks){
    
    vector<double> ans(2);
    double best_lnl = -DBL_MAX;

    for(int i = 0; i < x_ticks; i++){
        double site = (2.0 * i / (x_ticks-1) - 1) * abs(2.0 * i / (x_ticks-1) - 1) * width + center_point[0];
        for(int j = 0; j < y_ticks; j++){
            double f1 = (2.0 * j / (y_ticks-1) - 1) * abs(2.0 * j / (y_ticks-1) - 1) * height + center_point[1];
            
            vector<double> point(2);
            point[0] = site;
            point[1] = f1;
                
            double this_lnl = to_be_optimized_additive_fast(point);

            if(this_lnl > best_lnl) {
                best_lnl = this_lnl;
                ans = point;
            }
        }
    }

    return ans;
}

vector<double> selection_opt::grid_search_dominant0(vector<double> center_point, double width, double height, int x_ticks, int y_ticks){
    
    vector<double> ans(2);
    double best_lnl = -DBL_MAX;

    for(int i = 0; i < x_ticks; i++){
        double site = (2.0 * i / (x_ticks-1) - 1) * abs(2.0 * i / (x_ticks-1) - 1) * width + center_point[0];
        for(int j = 0; j < y_ticks; j++){
            double f1 = (2.0 * j / (y_ticks-1) - 1) * abs(2.0 * j / (y_ticks-1) - 1) * height + center_point[1];
            
            vector<double> point(2);
            point[0] = site;
            point[1] = f1;
                
            double this_lnl = to_be_optimized_pop0_dominant_fast(point);

            if(this_lnl > best_lnl) {
                best_lnl = this_lnl;
                ans = point;
            }
        }
    }

    return ans;
}

vector<double> selection_opt::grid_search_dominant1(vector<double> center_point, double width, double height, int x_ticks, int y_ticks){
    
    vector<double> ans(2);
    double best_lnl = -DBL_MAX;

    for(int i = 0; i < x_ticks; i++){
        double site = (2.0 * i / (x_ticks-1) - 1) * abs(2.0 * i / (x_ticks-1) - 1) * width + center_point[0];
        for(int j = 0; j < y_ticks; j++){
            double f1 = (2.0 * j / (y_ticks-1) - 1) * abs(2.0 * j / (y_ticks-1) - 1) * height + center_point[1];
            
            vector<double> point(2);
            point[0] = site;
            point[1] = f1;
                
            double this_lnl = to_be_optimized_pop1_dominant_fast(point);

            if(this_lnl > best_lnl) {
                best_lnl = this_lnl;
                ans = point;
            }
        }
    }

    return ans;
}

#endif