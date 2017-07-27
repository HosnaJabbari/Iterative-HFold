
// a simple driver for the pairfold

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

// include the pairfold library
#include "pairfold.h"
#include "externs.h"


int main (int argc, char *argv[])
{
    char sequence1[MAXSLEN], sequence2[MAXSLEN];
    char structure[2*MAXSLEN];
    double energy;
    int i;
    char structures[MAXSUBSTR][MAXSLEN];
    double energies[MAXSUBSTR];

    if (argc != 2)
    {
        printf ("Usage: %s <input_file>\n", argv[0]);
        printf ("Example: %s owczarzy04.txt\n");
        return 0;
    }

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
    int dna_or_rna = DNA;

    // temperature: any integer or real number between 0 and 100
    // represents degrees Celsius
    double temperature = 37;

    // initialize the thermodynamic parameters
    // call init_data only once for the same dna_or_rna and same temperature
    // if one of them changes, call init_data again
    init_data (argv[0], config_file, dna_or_rna, temperature);


    char buffer[256]; // read each line from file in buffer
    FILE *file;
    char tmpstr[MAXSLEN];
    double conc;
    conc = 2e-6;
    double Tm1, Tm2;

    // open file and check if operation was successful
    if ( (file = fopen (argv[1],"r")) == NULL )
    {
            fprintf(stderr, "Cannot open file: %s\n", argv[2]);
            exit(0);
    }
    // if file was opened, read first line
    fgets(buffer, sizeof(buffer), file);

    // read file line by line until you reach end of file
    while (!feof(file))
    {
        // skip comments and blank lines
        if (buffer[0] == '#' || buffer[0] == '\n')
        {
            fgets (buffer, sizeof(buffer), file);
            continue;
        }
        sscanf (buffer, "%s", tmpstr);

        Tm1 = calc_Tm_complementary (tmpstr, conc, conc, 1.0);
        Tm2 = calc_Tm_complementary_with_entropy_santalucia_2004 (tmpstr, conc, conc, 1.0);
        printf ("%30s\t%.2lf\t%.2lf\n", tmpstr, Tm1, Tm2);
        fgets (buffer, sizeof(buffer), file);
    }
    fclose (file);

    return 0;
}



