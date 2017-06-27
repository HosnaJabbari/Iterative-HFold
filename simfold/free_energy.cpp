
// a simple driver for the simfold

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

// include the simfold header files
#include "simfold.h"
#include "externs.h"
#include "constants.h"
#include "params.h"

int main (int argc, char *argv[])
{
    char sequence[MAXSLEN];
    char structure[MAXSLEN];
    char restricted[MAXSLEN];
    char paramfile[200];
    double energy;
    
    if (argc != 4 && argc != 5)
    {
        printf ("Usage: %s <sequence> <structure> <param_set> [<restricted>]\n", argv[0]);
        return 0;
    }        
    strcpy (sequence, argv[1]);    
    strcpy (structure, argv[2]);
    strcpy (paramfile, argv[3]);
    if (argc == 5)
    {
        strcpy (restricted, argv[4]);
    }
    
    // Before calling any function in the library, you have to initialize config_file, dna_or_rna, temperature
    //     and to call the function init_data, which loads the thermodynamic parameters into memory
    
    // configuration file, the path should relative to the location of this executable
    char config_file[200] = "params/1multirnafold.conf";

    // what to fold: RNA or DNA
    int dna_or_rna = RNA;

    // temperature: any integer or real number between 0 and 100
    // represents degrees Celsius
    double temperature = 37; 
    
    // initialize the thermodynamic parameters
    // call init_data only once for the same dna_or_rna and same temperature
    // if one of them changes, call init_data again   
    init_data (argv[0], config_file, dna_or_rna, temperature);         

    num_params = create_string_params ();
        
    // Before calling any function in the library, you have to initialize config_file, dna_or_rna, temperature
    //     and to call the function init_data, which loads the thermodynamic parameters into memory        
    fill_data_structures_with_new_parameters (paramfile);
    
    if (argc == 3)
    {
        energy = free_energy_simfold (sequence, structure);
    }    
    else
    {
        energy = free_energy_simfold_restricted (sequence, structure, restricted);
    }    
    printf ("%.2lf\n", energy);
    return 0;
}
    


