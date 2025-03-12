
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

//kevin July 17 2017
#include "hfold_validation.h"
#include <unistd.h>

void printUsage();

int main (int argc, char *argv[])
{
    char sequence[MAXSLEN];
    char structure[MAXSLEN];
    char restricted[MAXSLEN];
    double energy;
    char structures[MAXSUBSTR][MAXSLEN];
    double energies[MAXSUBSTR];


    //kevin: june 22 2017
	//validation for command line argument
    char* inputPath;
	inputPath = (char*) malloc(sizeof(char) * 1000);

	char* outputPath;
	outputPath = (char*) malloc(sizeof(char) * 1000);

	bool sequenceFound = false;
	bool restrictedFound = false;
	bool inputPathFound = false;
	bool outputPathFound = false;
	bool errorFound = false;
	int option;
	//kevin: june 22 2017 https://www.gnu.org/software/libc/manual/html_node/Example-of-Getopt.html#Example-of-Getopt
	while ((option = getopt (argc, argv, "s:r:i:o:")) != -1){
		switch (option)
		{
		case 's':
			if(sequenceFound){
				fprintf(stderr,"-s is duplicated\n");
				errorFound = true;
				break;
			}
			if(inputPathFound){
				fprintf(stderr,"Cannot combine -i with -s/-r \n");
				errorFound = true;
				break;
			}
			strcpy(sequence, optarg);
			//printf("seq: %s\n",sequence);
			sequenceFound = true;
			break;
		case 'r':
			if(restrictedFound || inputPathFound){
				fprintf(stderr,"-r is duplicated\n");
				errorFound = true;
				break;
			}
			if(inputPathFound){
				fprintf(stderr,"Cannot combine -i with -s/-r \n");
				errorFound = true;
				break;
			}
			strcpy(restricted, optarg);
			//printf("restricted: %s\n",restricted);
			restrictedFound = true;
			break;
		case 'i':
			if(restrictedFound || sequenceFound){
				fprintf(stderr,"Cannot combine -i with -s/-r \n");
				errorFound = true;
				break;
			}
			strcpy(inputPath,optarg);
			//printf("file: %s %d\n", file,access(file, F_OK|R_OK));
			if(access(inputPath, F_OK) == -1) { //if file does not exist
				fprintf(stderr,"Input file not exist\n");
				exit(4);
			}
			if (!validateHFOLDInputFile(inputPath, sequence, restricted)) {
				fprintf(stderr,"Input file is invalid\n");
				errorFound = true;
				break;
			}
			inputPathFound = true;
			break;
		case 'o':
			strcpy(outputPath, optarg);
			//printf("access: %d\n",access(outputPath, F_OK));
			if(access(outputPath, F_OK) != -1) { //if file already exist
				addTimestamp(&outputPath);
			}
			outputPathFound = true;
			break;
		default:
			errorFound = true;
			break;
		}
		//clean up when error
		if(errorFound){
            free(inputPath);
			free(outputPath);
			printUsage();
			exit(1);
		}
	}

	if(!inputPathFound){
		//if sequence or restricted is missing when input file is not present
		if(!(sequenceFound && restrictedFound)){
			fprintf(stderr, "-s/-r is missing\n");
			printUsage();
			exit(1);
		}
	}

	if(!validateSequence(sequence)){
		fprintf(stderr, "-s is invalid\n");
		//printUsage();
		exit(1);
	}

	if(!validateStructure(restricted, sequence)){
		fprintf(stderr, "-r is invalid\n");
		//printUsage();
		exit(1);
	}else{
		replaceBrackets(restricted);
	}

	//kevin: june 22 2017 if we have output path and input path, try to combine both
	if(outputPathFound && inputPathFound){
		addPath(&outputPath, inputPath);
		//printf("out path: %s\n",outputPath);
	}
	//kevin: june 22 2017
	//end of validation for command line arguments

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
            fprintf (stderr, "There is something wrong with the structure, doesn't match restricted\n");
			fprintf (stderr, "  %s\n  %s\n  %s\t%.2lf\n", sequence, restricted, structure, energy);
			exit(1);
        }
    }

    printf ("Seq: %s\n", sequence);
    printf ("RES: %s  %.2lf\n", structure, energy);

    return 0;
}



