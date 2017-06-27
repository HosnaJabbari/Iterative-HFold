
// a simple driver for the simfold

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

// include the simfold header files
#include "simfold.h"
#include "externs.h"
#include "constants.h"
#include "params.h"
#include "s_partition_function.h"

char sequence[MAXSLEN];
char parameter_filename[500] = "";
double threshold = 0.0;

void usage (const char *prog)
// PRE:  None
// POST: Prints the usage message and exits
{
    printf ("Usage: %s -s <sequence> [options]\n\n", prog);
    printf ("  -s <sequence>\n");
    printf ("\tThe RNA sequence to fold, max length is %d\n\n", MAXSLEN);
    printf ("Options:\n");
    printf ("  -t <probability_threshold>\n");
    printf ("\tDisplays probabilities above given threshold. DEFAULT %.2lf\n\n", threshold);
    printf ("  -p <parameter_file>\n");
    printf ("\tIf none specified, it is using the Turner99 parameters.\n\n");
    printf ("  -d\n\tTurn debug on, by default it is off\n");
    printf ("  -h\n\tPrint this help message\n\n");
    printf ("Examples:\n");
    printf ("\t%s -s \"GCAACGAUGACAUACAUCGCUAGUCGACGC\"\n", prog);
    printf ("\t%s -s \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" -t 0.5\n", prog);
    printf ("\t%s -s \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" -p params/CG_best_parameters_ISMB2007.txt\n");
    exit (0);
}


void get_arguments (int argc, char *argv[])
// PRE:  None
// POST: Get the parameters from the command line
{
    int c;
    extern char *optarg;
    extern int optind;
    int errflag = 0;

    while ((c = getopt (argc, argv, "s:t:p:dh?")) != -1)
    {
        switch (c)
        {
            case 's':
                strcpy (sequence, optarg);
                break;
            case 't':
                threshold = atof (optarg);
                break;                
            case 'p':
                strcpy (parameter_filename, optarg);
                break;
            case 'd':
                debug = 1;
                break;
            case 'h':
            case '?':
            default:
                errflag = 1;
        }
    }
    if (errflag || optind != argc || argc < 3)
        usage (argv[0]);

}



int main (int argc, char *argv[])
{
    char structure[MAXSLEN];    
    double energy;
    char structures[MAXSUBSTR][MAXSLEN];
    double energies[MAXSUBSTR];    
    int actual_num_str;
    //debug = 1;

    get_arguments (argc, argv);
    
    // Before calling any function in the library, you have to initialize config_file, dna_or_rna, temperature
    //     and to call the function init_data, which loads the thermodynamic parameters into memory
    
    // configuration file, the path should be relative to the location of this executable
    char config_file[200];
    //strcpy (config_file, "params/multirnafold.conf");

    //kevin 26 june 2017
    //occurs in /simfold/melting_temp.cpp, /simfold/simfold.cpp, /simfold/simfold_pf.cpp, /simfold/test_get_counts.cpp
    //changed strcpy to str cat beacause with installed simfold, it is trying to look at /simfold/.libs/params/multirnafold.conf
    //PARAMS_BASE_PATH is in constants.h to be ../params/ to get away from .libs directory
    strcat(config_file, PARAMS_BASE_PATH);
    strcat(config_file, "multirnafold.conf");

    // what to fold: RNA or DNA
    int dna_or_rna;
    dna_or_rna = RNA;

    // temperature: any integer or real number between 0 and 100
    // represents degrees Celsius
    double temperature = 37.0; 
    
    // initialize the thermodynamic parameters
    // call init_data only once for the same dna_or_rna and same temperature
    // if one of them changes, call init_data again   
    init_data (argv[0], config_file, dna_or_rna, temperature);         

    if (strlen (parameter_filename) > 0)
        fill_data_structures_with_new_parameters (parameter_filename);
    else
        fill_data_structures_with_new_parameters ("params/turner_parameters_fm363_constrdangles.txt");

    printf ("Seq: %s\n", sequence);

    // create the partition function object. By default, it includes dangling ends.
    //  Add 1 to the constructor if you want not to include the dangling ends (faster)
    // WARNING: if your sequence is longer than 1000, you may get overflow. This will be fixed in a next release.
    s_partition_function *part = new s_partition_function (sequence);

    // compute the partition function
    double pf = part->compute_partition_function();
    printf ("Partition function = %g\n", pf);

    // compute the base pair probabilities if you need (the partition function has to be computed)
    part->compute_base_pair_probabilities();

    // save the dot plot file (the lower triangle is empty for now)
    part->PS_dot_plot("dot.ps");


    // You can print the base pair probabilities above some threshold, for example 0.9
    part->print_base_pair_probabilities(threshold);

    delete part;
    return 0;
}
    
  
