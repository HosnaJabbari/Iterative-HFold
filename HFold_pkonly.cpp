
// a simple driver for the HFold

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

// include the simfold header files
#include "simfold.h"
#include "externs.h"
#include "h_globals.h"
//#include "h_externs.h"
#include "constants.h"
#include "params.h"


// Hosna June 20th, 2007
//#include "W_final.h"
#include "hfold_pkonly.h"

int main (int argc, char *argv[])
{
    char sequence[MAXSLEN];
    char structure[MAXSLEN];
    char restricted[MAXSLEN];
    double energy;
    char structures[MAXSUBSTR][MAXSLEN];
    double energies[MAXSUBSTR];

    if (argc != 3)
    {
        printf ("Usage: %s <sequence> <restricted_structure>\n", argv[0]);
        printf ("Example: %s \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" \"(____________________________)\"\n", argv[0]);
        printf ("\tRestricted structure symbols:\n");
        printf ("\t\t() restricted base pair\n");
        printf ("\t\t_ no restriction\n");
        return 0;
    }
    strcpy (sequence, argv[1]);
    strcpy (restricted, argv[2]);

    // Before calling any function in the library, you have to initialize config_file, dna_or_rna, temperature
    //     and to call the function init_data, which loads the thermodynamic parameters into memory

    // configuration file, the path should be relative to the location of this executable
    char config_file[200];
    strcpy (config_file, "./simfold/params/multirnafold.conf");

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

	// Hosna, July 18, 2012
	// In simfold we have the following for RNA && temp=37
	fill_data_structures_with_new_parameters ("./simfold/params/turner_parameters_fm363_constrdangles.txt");

	// Hosna, July 25, 2012
	// in HotKnots and ComputeEnergy package the most up-to-date parameters set is DP09.txt
	// so we add it here
	fill_data_structures_with_new_parameters ("./simfold/params/parameters_DP09.txt");
	
	energy = hfold_pkonly(sequence, restricted, structure); 

    //delete min_fold;

    // check if restricted is included in structure
	// Hosna March 7, 2012
	// for optimality we get the value once, use it many times
	int seqLen = strlen (sequence);
    for (int i=0; i < seqLen; i++)
    {
        if ((restricted[i] == '(' || restricted[i] == ')' || restricted[i] == '.') &&
            (restricted[i] != structure[i]))
        {
            printf ("There is something wrong with the structure, doesn't match restricted\n");
			printf ("  %s\n  %s\n  %s\t%.2lf\n", sequence, restricted, structure, energy);
			exit(1);
        }
    }

    printf ("Seq: %s\n", sequence);
    printf ("RES: %s  %.2lf\n", structure, energy);

    return 0;
}



