
// test the function get_feature_counts_restricted

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

// include the simfold header files
#include "simfold.h"
#include "externs.h"
#include "constants.h"
#include "params.h"

/*
int check_counts_linear (int numpars, double *params, double *c, double f, double energy)
// Return 1 if energy ~= c'x + f    (where x is params)
// Return 0 otherwise
// Added on Feb 29, 2008
{
    double energy_c = 0;
    int i;
    for (i=0; i < numpars; i++)
    {
        energy_c +=  c[i] * params[i];  
    }
    energy_c += f;
    
    if (fabs(energy_c - energy) > 0.1)  // not really sure what's a good threshold. Maybe 0 is fine
    {
        printf ("ERROR! Something is wrong with the counts or the free energy: c'x+f = %.2lf, energy = %.2lf diff=%.2lf\n", energy_c, energy, fabs(energy_c-energy));
        return 0;
    }        
    return 1;
}
*/

int main (int argc, char *argv[])
{
    char sequence[MAXSLEN];
    char structure[MAXSLEN];
    char restricted[MAXSLEN];
    char paramfile[200];
    double energy;
    int ignore_dangles;
    int with_counts;
    
    if (argc != 5)
    {
        printf ("Usage: %s <sequence> <structure> <ignore_dangles:0|1> <with_counts:0|1>\n", argv[0]);
        printf ("Example: %s \"ACCAACAAAACAGG\" \".((..<xxxx>.))\" 0 1\n", argv[0]); 
        return 0;
    }        
    strcpy (sequence, argv[1]);    
    strcpy (structure, argv[2]);
    ignore_dangles = atoi (argv[3]);
    with_counts = atoi (argv[4]);
    debug = 0;
    
    // Before calling any function in the library, you have to initialize config_file, dna_or_rna, temperature
    //     and to call the function init_data, which loads the thermodynamic parameters into memory
    
    // configuration file, the path should relative to the location of this executable
    //char config_file[200] = "params/multirnafold.conf";
    
    //kevin 26 june 2017
    //occurs in /simfold/melting_temp.cpp, /simfold/simfold.cpp, /simfold/simfold_pf.cpp, /simfold/test_get_counts.cpp
    //changed strcpy to str cat beacause with installed simfold, it is trying to look at /simfold/.libs/params/multirnafold.conf
    //PARAMS_BASE_PATH is in constants.h to be ../params/ to get away from .libs directory
    char config_file[200];
    strcat(config_file, PARAMS_BASE_PATH);
    strcat(config_file, "multirnafold.conf");

    // what to fold: RNA or DNA
    int dna_or_rna = RNA;

    // temperature: any integer or real number between 0 and 100
    // represents degrees Celsius
    double temperature = 37; 
    
    // initialize the thermodynamic parameters
    // call init_data only once for the same dna_or_rna and same temperature
    // if one of them changes, call init_data again   
    init_data (argv[0], config_file, dna_or_rna, temperature);

    // I have to call the following function because the counter function doesn't do minimization of dangling ends
    //fill_data_structures_with_new_parameters ("params/turner_parameters_fm363_constrdangles.txt");
    fill_data_structures_with_new_parameters ("params_sub-1.txt");
    
    int num_parameters = get_num_params();
    PARAMTYPE params[num_parameters];
    double dparams[num_parameters];
    double c[num_parameters];
    double f;

    if (with_counts)
    {
        save_parameters_in_array (params);
        // this function saves the params*100, so divide by 100 here
        for (int i=0; i < num_parameters; i++)  dparams[i] = params[i]/100.0;
        energy = get_feature_counts_restricted (sequence, structure, c, f, 1, ignore_dangles, 1);
        //energy = get_feature_counts (sequence, structure, c, f);
        print_counter (c, f);
        if (check_counts_linear (num_parameters, dparams, c, f, energy))
        {
            printf ("Energy and c'x+f MATCH!\n");
        }
    }
    else
    {
        energy = get_feature_counts_restricted (sequence, structure, NULL, f, 1, ignore_dangles, 0);        
    }
    printf ("Energy is %.2lf\n", energy);
        
    return 0;
}
    


