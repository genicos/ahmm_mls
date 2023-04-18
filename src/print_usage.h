#ifndef __PRINT_USAGE_H
#define __PRINT_USAGE_H

void print_usage() {
    
    cerr << endl << endl << "ahmm_mls usage:" << endl << endl ;
    cerr << "\trequired:" << endl ;
    cerr << "\t\t-i [string]\t\tinput file" << endl ;
    cerr << "\t\t-s [string]\t\tsample id and ploidy file" << endl ;
    cerr << "\t\t-m [string]\t\tmodel file" << endl ;

    cerr << "\toptional:" << endl ;
    cerr << "\t\t--help\t\t\tprint this help statement" << endl ;
    cerr << "\t\t-g\t\t\tsamples are specified with genotypes rather than read counts" << endl ;
    cerr << "\t\t-c [int]\t\tnumber of cores used" << endl ;
    cerr << "\t\t-f\t\t\toutput relative fitnesses rather than selection coefficients" << endl ;
    cerr << "\t\t-vo\t\t\tverbose stderr output" << endl ;
    cerr << "\t\t-R [float]\t\tspecify morgan distance from selected site after which they are ignored" << endl ;
    cerr << "\t\t-t1 [float]\t\tthreshold of lnL ratio range in simplex for first stage of optimization" << endl ;
    cerr << "\t\t-t2 [float]\t\tthreshold of lnL ratio range in simplex for second stage of optimization" << endl ;
    cerr << "\t\t-k [int]\t\tnumber of skipped regions between adjacent sampled sites\n\t\t\t\t\tfor each one calculated." << endl ;
}

#endif

