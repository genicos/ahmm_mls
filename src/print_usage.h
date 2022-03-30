#ifndef __PRINT_USAGE_H
#define __PRINT_USAGE_H

void print_usage() {
    
    cerr << endl << endl << "ahmm_mls usage:" << endl << endl ;
    cerr << "\trequired:" << endl ;
    cerr << "\t\t-i [string]\t\tinput file name" << endl ;
    cerr << "\t\t-s [string]\t\tsample id and ploidy file" << endl ;
    cerr << "\t\t-m [float] [int]\t\tadmixture information with admixture fraction, time" << endl ;

    cerr << "\tat least one required:" << endl;
    cerr << "\t\t-l [string]\t\tsite file" << endl ;
    cerr << "\t\t-M [string]\t\tmodel file" << endl ;

    cerr << "\toptional:" << endl ;
    cerr << "\t\t--help\t\t\tprint this help statement" << endl ;
}

#endif

