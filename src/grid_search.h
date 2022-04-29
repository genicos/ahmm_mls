#ifndef GRID_SEARCH
#define GRID_SEARCH
#include <vector>
#include <math.h>

vector<double> selection_opt::grid_search(vector<double> center_point, double width, double height, double x_point_sep, double y_point_sep){
    vector<double> ans(3);
    double best_lnl = -MAX_DBL;

    for(double site = center_point[0] - width/2; site >= 0 && site <= center_point[0] + width/2 && site <= chromosome_size; site += x_point_sep){
        for(double f1 = center_point[1] - height/2; f1 >= 0 && f1 <= center_point[1] + height/2; f1 += y_point_sep){
            for(double f2 = center_point[2] - height/2; f2 >= 0 && f2 <= center_point[2] + height/2; f2 += y_point_sep){
                vector<double> point(3);
                point[0] = site;
                point[1] = f1;
                point[2] = f2;
                
                double this_lnl = to_be_optimized_only_near_sites(point);
                if(this_lnl > best_lnl){
                    best_lnl = this_lnl
                    ans = point;
                }
            }
        }
    }

    return ans;
}

#endif