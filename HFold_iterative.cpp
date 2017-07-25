// a simple driver for the HFold
#include <iostream>
#include <algorithm>
#include <iterator>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <stack>
#include <ctime>
#include <dirent.h>
#include <signal.h>
#include <pthread.h>

// include the simfold header files
#include "simfold.h"
#include "externs.h"
#include "h_globals.h"
//#include "h_externs.h"
#include "constants.h"
#include "params.h"
#include "common.h"

#include <sys/types.h>
#include <sys/wait.h>
#include <string>

// Hosna June 20th, 2007
//#include "W_final.h"
//#include "hfold.h"
#include "hfold_iterative.h"

#include <ctime>
#include <regex>
#include <regex.h>

#include "hfold_validation.h" //kevin June 22 2017

#define HFOLD 						"./HFold"
#define HFOLD_PKONLY 				"./HFold_pkonly"
#define HFOLD_INTERACTING 			"./HFold_interacting"
#define HFOLD_INTERACTING_PKONLY 	"./HFold_interacting_pkonly"
#define SIMFOLD 					"./simfold"

#define INPUTPATH 			"./it_test/exons/exon_6_mutated_flank_input_PMO_1/"
#define LOGFILEPATH 		"./logfile.txt"
#define OUTPUTPATH 			"./it_test/exons/exon_6_mutated_flank_output_PMO/"
#define OUTPUTPATH_HFOLD	"./temp/"


FILE *logFile;
char *file;
int timeout = 1200;
static std::smatch match;
static std::regex sequenceRegex;
static std::regex structureRegex;
static std::regex replacedRegex;
static std::string sequenceString;
static std::string structureString;
static std::string replacementText;
static std::string energyString;
static std::stack<int> customStack;


/* I am changing iterative HFold so that we run 4 methods and choose the structure with the lowest energy as the winner
1) run HFold_PKonly on input structure, if pked, keep the pk bases as input structure and run HFold on them and get the result
2) run HFold on the original input and get the result
3) get the beginning and end of the given structure, give the subsequence and the substructure to simfold restricted as input, get the simfold structure, then give that as input to HFold_PKonly, get the pked bases as new input to give HFold and get the result
4) run simfold restricted with the given input sequence and structure, and get the simfold structure. Only choose part of the structure that contains the original given input structure, then give that to HFold_PKonly as input, get pked bases and run HFold on them and get the result. */
int main (int argc, char **argv) {

	void *res;

	char sequence[MAXSLEN];
	char structure[MAXSLEN];

	char *output_path;
	char *method1_structure = (char*) malloc(sizeof(char) * MAXSLEN);
	char *method2_structure = (char*) malloc(sizeof(char) * MAXSLEN);
	char *method3_structure = (char*) malloc(sizeof(char) * MAXSLEN);
	char *method4_structure = (char*) malloc(sizeof(char) * MAXSLEN);
	char final_structure[MAXSLEN];

	double *method1_energy = (double*) malloc(sizeof(double) * INF);
	double *method2_energy = (double*) malloc(sizeof(double) * INF);
	double *method3_energy = (double*) malloc(sizeof(double) * INF);
	double *method4_energy = (double*) malloc(sizeof(double) * INF);
	double final_energy = INF;

	int files_length = -1;
        int method_chosen = -1;
	 
	int **result;
	char **files_found;
	
        result  = (int**) malloc(sizeof(int*) * MAXSLEN);
        for (int i=0; i < MAXSLEN; i++) {
                result[i] = (int*) malloc(sizeof(int) * 2);
        }



	output_path = (char*) malloc(sizeof(char) * 1000);
	file = (char*) malloc(sizeof(char) * 1000);
	
	logFile = fopen(LOGFILEPATH, "w");
	if (logFile == NULL) {
		giveup("Could not open logfile", LOGFILEPATH);
	}

	//kevin: june 22 2017
	//validation for command line argument
	bool sequenceFound = false;
	bool structureFound = false;
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
				printf("-s is duplicated\n");
				errorFound = true;
				break;
			}
			if(inputPathFound){
				printf("Cannot combine -i with -s/-r \n");
				errorFound = true;
				break;
			}
			strcpy(sequence, optarg);
			//printf("seq: %s\n",sequence);
			sequenceFound = true;
			break;
		case 'r':
			if(structureFound || inputPathFound){
				printf("-r is duplicated\n");
				errorFound = true;
				break;
			}
			if(inputPathFound){
				printf("Cannot combine -i with -s/-r \n");
				errorFound = true;
				break;
			}
			strcpy(structure, optarg);
			//printf("structure: %s\n",structure);
			structureFound = true;
			break;
		case 'i':
			if(structureFound || sequenceFound){
				printf("Cannot combine -i with -s/-r \n");
				errorFound = true;
				break;
			}
			strcpy(file,optarg);
			//printf("file: %s %d\n", file,access(file, F_OK|R_OK));
			if(access(file, F_OK) == -1) { //if file does not exist
				printf("Input file not exist\n");
				exit(4);
			}	
			if (!get_sequence_structure(file, sequence, structure)) {
				printf("Input file is invalid\n");
				errorFound = true;
				break;
			}
			inputPathFound = true;
			break;
		case 'o':
			strcpy(output_path, optarg);
			//printf("access: %d\n",access(output_path, F_OK));
			if(access(output_path, F_OK) != -1) { //if file already exist
				addTimestamp(&output_path);
			}	
			outputPathFound = true;
			break;
		default:
			errorFound = true;
			break;
		}
		//clean up when error
		if(errorFound){
			printUsage();
			exit(1);
		}
	}

	if(!inputPathFound){
		//if sequence or structure is missing when input file is not present
		if(!(sequenceFound && structureFound)){
			printf("-s/-r is missing\n");
			printUsage();
			exit(1);
		}
	}

	if(!validateSequence(sequence)){
		printf("-s is invalid\n");
		printUsage();
		exit(1);
	}

	if(!validateStructure(structure, sequence)){
		printf("-r is invalid\n");
		printUsage();
		exit(1);
	}else{
		replaceBrackets(structure);
	}

	//if we have output path and input path, try to combine both
	if(outputPathFound && inputPathFound){
		addPath(&output_path, file);
		//printf("out path: %s\n",output_path);
		
	}	
	//kevin: june 22 2017
	//end of validation for command line arguments

	write_log_file("Starting Program", "", 'I');


	*method1_energy = INF;
	*method2_energy = INF;
	*method3_energy = INF;
	*method4_energy = INF;

	method1_structure[0] = NULL;
	method2_structure[0] = NULL;
	method3_structure[0] = NULL;
	method4_structure[0] = NULL;
	final_structure[0] = NULL;		


        /*************** First Method ***************/
        method1_calculation(sequence, structure, method1_structure, method1_energy);

        /*************** Second Method ***************/
        method2_calculation(sequence, structure, method2_structure, method2_energy);

        /*************** Third Method ***************/
        method3_calculation(sequence, structure, method3_structure, method3_energy);

        /*************** Fourth Method ***************/
        method4_calculation(sequence, structure, method4_structure, method4_energy, result);


        if (*method1_energy < final_energy && *method1_energy != 0) {
                final_energy = *method1_energy;
                strcpy(final_structure, method1_structure);
                method_chosen = 1;
        }

        if (*method2_energy < final_energy && *method2_energy != 0) {
                final_energy = *method2_energy; 
                strcpy(final_structure, method2_structure);             
                method_chosen = 2;
        }

        if (*method3_energy < final_energy && *method3_energy != 0) {
                final_energy = *method3_energy; 
                strcpy(final_structure, method3_structure);             
                method_chosen = 3;
        }

        if (*method4_energy < final_energy && *method4_energy != 0) {
                final_energy = *method4_energy; 
                strcpy(final_structure, method4_structure);             
                method_chosen = 4;
        }


	if (final_energy == INF || method_chosen == -1) {
		write_log_file("Could not find energy", "", 'E');
		printf("ERROR: could not find energy\n");
		final_energy = NULL;
		exit(6);
	}
	

	
	
	//kevin: june 22 2017
	//output to file
	if(outputPathFound){
		if(!save_file("", output_path, sequence, structure, final_structure, final_energy, method_chosen)){
			printf("write to file fail\n");
			exit(4);
		}
	}else{
		//kevin: june 22 2017
		//changed format for ouptut to stdout
		std::cout << "Seq: " << sequence << "\n";
		std::cout << "RES: " << final_structure << "  " << final_energy << "\n" << std::flush;
		//std::cout << "\nfinal structure: " << final_structure << " final energy: " << final_energy << "\n\n" << std::flush;
	}

	

	write_log_file("Ending Program", "", 'I');
	
	// Clean up

	free(file);
           	
	//free(tinfo);
	free(method1_structure);
	free(method2_structure);
	free(method3_structure);
	free(method4_structure);
	free(method1_energy);
	free(method2_energy);
	free(method3_energy);
	free(method4_energy);
        free(output_path);
        for (int i=0; i < MAXSLEN; i++) {
                free(result[i]);
        }
	free(result);


	// This will cause a seg fault if you have the log file open and are watching it.
	if (logFile) {
		fclose(logFile);
	}

	return 0;
}



void printUsage(){
	/*
	printf ("\nUsage: HFold_iterative <sequence> <structure>\n");
	printf ("Example: ./HFold_iterative \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" \"(____________________________)\" \n");
	printf ("Or \nUsage: HFold_iterative <path to input file> <path to output file>\n");
	printf ("Example: ./HFold_iterative \"/home/username/Desktop/inputFile.txt\" -o \"/home/username/Desktop/outFile.txt\" \n");
	printf ("\tRestricted structure symbols:\n");
	printf ("\t\t() restricted base pair\n");
	printf ("\t\t _ no restriction\n");
*/
	printf("Usage ./HFold_iterative -s <sequence> -r <structure> [-o </path/to/file>]\n");
	printf("or\n");
	printf("Usage ./HFold_iterative -i </path/to/file> [-o </path/to/file>]\n");
	printf ("  Restricted structure symbols:\n");
	printf ("    () restricted base pair\n");
	printf ("    _ no restriction\n");
	printf("Example:\n");
	printf("./HFold_iterative -s \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" -r \"(____________________________)\"\n");
	printf("./HFold_iterative -i \"/home/username/Desktop/myinputfile.txt\" -o \"/home/username/Desktop/some_folder/outputfile.txt\"\n");
	printf("Please read README for more details\n");
}


void segfault_sigaction(int signal, siginfo_t *si, void *arg) {
	//Need to expand on this.
    printf("Caught segfault at address %p\n", si->si_addr);
    exit(si->si_errno);
}


void method1_calculation (char *sequence, char *structure, char *method1_structure, double *method1_energy) {
	char new_input_structure[MAXSLEN] = "\0";	
	char hfold_structure[MAXSLEN] = "\0";	
	char hfold_pkonly_structure[MAXSLEN] = "\0";	
	
	double hfold_energy = 0;	
	double hfold_pkonly_energy = 0;

	bool has_pk;

	if (!call_HFold(HFOLD_PKONLY, sequence, structure, hfold_pkonly_structure, &hfold_pkonly_energy)) {
		*method1_energy = hfold_pkonly_energy;		
		return;
	}

	has_pk = find_new_structure(hfold_pkonly_structure, new_input_structure);
	//std::cout << "hfold_pkonly_structure = " << hfold_pkonly_structure << " new_input_structure = " << new_input_structure << " has_pk = " << has_pk << '\n' << std::flush;

	if (has_pk) {		
		if (!call_HFold(HFOLD, sequence, new_input_structure, hfold_structure, &hfold_energy)) {
			std::cout << "HFold failed " << std::endl;
			*method1_energy = hfold_energy;
			return;
		}
		*method1_energy = hfold_energy;
		strcpy(method1_structure, hfold_structure);

	} else {
		*method1_energy = hfold_pkonly_energy;
		strcpy(method1_structure, hfold_pkonly_structure);
	}

	//std::cout << "method1_structure = " << method1_structure << " method1_energy = " << *method1_energy << '\n' << std::flush;

	if (method1_structure == "") {
		write_log_file("The structure should not be null.", file, 'E');
	}
}

void method2_calculation (char *sequence, char *structure, char *method2_structure, double *method2_energy) {
	char hfold_structure[MAXSLEN] = "\0";	
	
	double hfold_energy = 0;	
	
	if (!call_HFold(HFOLD, sequence, structure, hfold_structure, &hfold_energy)) {
		*method2_energy = hfold_energy;
		return;
	}
	*method2_energy = hfold_energy;
	strcpy(method2_structure, hfold_structure);

	//std::cout << "method2_structure = " << method2_structure << " method2_energy = " << *method2_energy << '\n' << std::flush;

	if (method2_structure == "") {
		write_log_file("The structure should not be null.", file, 'E');
	}
}


void method3_calculation (char *sequence, char *structure, char *method3_structure, double *method3_energy) {
	char sub_sequence[MAXSLEN] = "\0";
	char sub_structure[MAXSLEN] = "\0";
	
	char replaced_structure[MAXSLEN] = "\0";
	char new_input_structure[MAXSLEN] = "\0";
	char hfold_structure[MAXSLEN] = "\0";
	char hfold_pkonly_structure[MAXSLEN] = "\0";
	char simfold_structure[MAXSLEN] = "\0";

	double hfold_energy = 0;	
	double hfold_pkonly_energy = 0;
	double simfold_energy = 0;
	
	int begin = -1, end = 0;
	bool has_pk;

	find_sub_sequence_structure(sequence, structure, sub_sequence , sub_structure, &begin, &end);


//	std::cout << "found sub sequence = " << sub_sequence << " and sub structure = " << sub_structure << "  with begin = " << begin << " and end = " << end << '\n' << std::flush;

	if (!call_simfold(SIMFOLD, sub_sequence, sub_structure, simfold_structure, &simfold_energy)) {
//                printf("TEST, method3_calculation's calling simfold result: %lf\n", simfold_energy);
		*method3_energy = simfold_energy;		
		return;
	}

	replace_simfold_partial_structure_with_original(structure, simfold_structure, replaced_structure, begin, end);
	//std::cout << "the replaced string is: " << replaced_structure << '\n' << std::flush;

	if (!call_HFold(HFOLD_PKONLY, sequence, replaced_structure, hfold_pkonly_structure, &hfold_pkonly_energy)) {
		*method3_energy = hfold_pkonly_energy;
		return;
	}
	has_pk = find_new_structure(hfold_pkonly_structure, new_input_structure);
	//std::cout << "hfold_pkonly_structure = " << hfold_pkonly_structure << " new_input_structure = " << new_input_structure << " has_pk = " << has_pk << '\n' << std::flush;

	if (has_pk) {
		if (!call_HFold(HFOLD, sequence, new_input_structure, hfold_structure, &hfold_energy)) {
			*method3_energy = hfold_energy;
			return;
		}
		*method3_energy = hfold_energy;
		strcpy(method3_structure, hfold_structure);

	} else {
		*method3_energy = hfold_pkonly_energy;
		strcpy(method3_structure, hfold_pkonly_structure);
	}

	//std::cout << "method3_structure = " << method3_structure << " method3_energy = " << *method3_energy << '\n' << std::flush;

	if (method3_structure == "") {
		write_log_file("The structure should not be null.", file, 'E');
	}
}

void method4_calculation (char *sequence, char *structure, char *method4_structure, double *method4_energy, int **result) {
	// we know the beginning and end of the given input substructure in the structure
	// now we need to run simfold restricted with the given structure, find the stem which contains the given substructure
	// then only include this structure as input structure, and do as we did in method 3
	char sub_sequence[MAXSLEN] = "\0";
	char sub_structure[MAXSLEN] = "\0";
	
	char replaced_structure[MAXSLEN] = "\0";
	char new_input_structure[MAXSLEN] = "\0";
	char hfold_structure[MAXSLEN] = "\0";
	char hfold_pkonly_structure[MAXSLEN] = "\0";
	char simfold_structure[MAXSLEN] = "\0";

	double hfold_energy = 0;	
	double hfold_pkonly_energy = 0;
	double simfold_energy = 0;
	
	int begin = -1, end = 0, B = 0, Bp = 0;
	int result_length = 0;
	bool has_pk;	

	find_sub_sequence_structure(sequence, structure, sub_sequence , sub_structure, &begin, &end);	
	strcpy(replaced_structure, structure);
	
	if (!call_simfold(SIMFOLD, sequence, structure, simfold_structure, &simfold_energy)) {
		*method4_energy = simfold_energy;
		return;
	}

	if (!find_each_substructure(structure, begin, end, result)) {
		return;
	}

	result_length = sizeof(result) / sizeof(result[0]);
	strcpy(replaced_structure, structure);
	for (int i = 0; i < result_length; i++) {
		find_independant_structures (simfold_structure, result[i][0], result[i][1], &B, &Bp);
		//std::cout << "found B = " << B << " and Bp = " << Bp << '\n' << std::flush;
		replace_simfold_structure_with_original(replaced_structure, simfold_structure, B, Bp);
	}
	
	//std::cout << "the replaced structure is: " << replaced_structure << '\n' << std::flush;

	if (!call_HFold(HFOLD_PKONLY, sequence, replaced_structure, hfold_pkonly_structure, &hfold_pkonly_energy)) {
		*method4_energy = hfold_pkonly_energy;
		return;
	}
	has_pk = find_new_structure(hfold_pkonly_structure, new_input_structure);
	//std::cout << "hfold_pkonly_structure = " << hfold_pkonly_structure << " new_input_structure = " << new_input_structure << " has_pk = " << has_pk << '\n' << std::flush;

	if (has_pk) {
		if (!call_HFold(HFOLD, sequence, new_input_structure, hfold_structure, &hfold_energy)) {
			*method4_energy = hfold_energy;
			return;
		}

		*method4_energy = hfold_energy;
		strcpy(method4_structure, hfold_structure);

	} else {
		*method4_energy = hfold_pkonly_energy;
		strcpy(method4_structure, hfold_pkonly_structure);
	}

	//std::cout << "method4_structure = " << method4_structure << " method4_energy = " << *method4_energy << '\n' << std::flush;

	if (method4_structure == "") {
		write_log_file("The structure should not be null.", file, 'E');
	}
}

//Calling HFold: programPath = HFOLD
//Calling HFold_PKonly: programPath = HFOLD_PKONLY
bool call_HFold (char *programPath, char *input_sequence, char *input_structure, char *output_structure, double *output_energy) {
	char config_file[200];
	strcpy (config_file, "./simfold/params/multirnafold.conf");

	//what to fold: RNA or DNA
	int dna_or_rna;
	dna_or_rna = RNA;
	//temperature: any integer or real number between 0 and 100
	// represents degrees Celsius
	double temperature = 37.0;

	// initialize the thermodynamic parameters
	// call init_data only once for the same dna_or_rna and same temperature
	// if one of them changes, call init_data again

	//init_data ("./HFold", config_file, dna_or_rna, temperature);
	init_data (programPath, config_file, dna_or_rna, temperature);
	
	fill_data_structures_with_new_parameters ("./simfold/params/turner_parameters_fm363_constrdangles.txt");

	// in HotKnots and ComputeEnergy package the most up-to-date parameters set is DP09.txt
	// so we add it here
	fill_data_structures_with_new_parameters ("./simfold/params/parameters_DP09.txt"); 
        
        if (strcmp(programPath, HFOLD) == 0) {
            *output_energy = hfold(input_sequence, input_structure, output_structure);
        }
        else if (strcmp(programPath, HFOLD_PKONLY) == 0){
                *output_energy = hfold_pkonly(input_sequence, input_structure, output_structure); 
        }else{
                printf("Error: invalid arguments are given: %s \nValid aurgumnets are: HFOLD and HFOLD_PKONLY\n");
                return false;
         }
	return true;

}

bool call_simfold (char *programPath, char *input_sequence, char *input_structure, char *output_structure, double *output_energy) {
        std::string result = "";	
	
	char config_file[200] = "simfold/";
	strcat(config_file, PARAMS_BASE_PATH);
	strcat(config_file, "multirnafold.conf");
	double temperature;
	temperature = 37; 
	init_data ("./simfold", config_file, RNA, temperature);         

        fill_data_structures_with_new_parameters ("simfold/params/turner_parameters_fm363_constrdangles.txt");
	// when I fill the structures with DP09 parameters, I get a segmentation fault for 108 base sequence!!!!
	// So I chopped the parameter set to only hold the exact number as the turner_parameters_fm363_constrdangles.txt, 
	// but still getting seg fault!
	fill_data_structures_with_new_parameters ("simfold/params/parameters_DP09_chopped.txt");

	printf ("Seq: %s\n", input_sequence);

	*output_energy = simfold_restricted (input_sequence, input_structure, output_structure);
//	printf ("Call_Simfold_RES( can be called by different methods): %s  %.2lf\n", output_structure, output_energy);
	return true;
}

void replace_simfold_partial_structure_with_original (char *input_structure, char *simfold_structure, char *replaced_structure, int begin, int end) {
	replacementText = "_";
	replacedRegex = "\\.";

	strcpy(replaced_structure, input_structure);
	int j = 0;
	
	for (int i = begin; i <= end; i++) {
		replaced_structure[i] = simfold_structure[j++];
	}

	structureString.assign(replaced_structure, strlen(replaced_structure));
	std::string regexResult(std::regex_replace(structureString, replacedRegex, replacementText));
	strcpy(replaced_structure, regexResult.c_str());
}

void replace_simfold_structure_with_original (char *replaced_structure, char *simfold_structure, int begin, int end) {
	replacementText = "_";
	replacedRegex = "\\.";

	for (int i = begin; i <= end; i++) {
		replaced_structure[i] = simfold_structure[i];
	}

	structureString.assign(replaced_structure, strlen(replaced_structure));
	std::string regexResult(std::regex_replace(structureString, replacedRegex, replacementText));
	strcpy(replaced_structure, regexResult.c_str());
}

bool get_sequence_structure (char *fileName, char *sequence, char *structure) {
	FILE *ioFile;
	size_t len = 0;
	ssize_t read = 0;
	char *sequenceBuffer = NULL;
	char *structureBuffer = NULL;

	replacementText = "";
	sequenceRegex = "(\\s|\n|\r)";
	structureRegex = "(\\s|\n|\r|\Z)";

	// Get structure and sequence from file
	ioFile = fopen(fileName, "r");
	if (ioFile != NULL) {
		if((read = getline(&sequenceBuffer, &len, ioFile)) == -1) {
			write_log_file("Could not read the first sequence", fileName, 'E');
			return false;
		}
		sequenceString.assign(sequenceBuffer, read);		

		if((read = getline(&structureBuffer, &len, ioFile)) == -1) {
			write_log_file("Could not read the first structure", fileName, 'E');
			return false;
		}
		structureString.assign(structureBuffer, read);

	} else {
		write_log_file("Could not open file", fileName, 'E');
		return false;
	}

	// Format sequence
	sequenceString = std::regex_replace(sequenceString, sequenceRegex, replacementText);
	std::cout << "Sequence:  " << sequenceString << " Length: " << sequenceString.length() << '\n' << std::flush;

	// Format structure
	structureString = std::regex_replace(structureString, structureRegex, replacementText);
	std::cout << "Structure: " << structureString << " Length: " << structureString.length() << '\n' << std::flush;

	// Copy strings
	strcpy(sequence, sequenceString.c_str());
	strcpy(structure, structureString.c_str());

	// Clean up
	fclose(ioFile);
	if (sequenceBuffer)
        free(sequenceBuffer);
	if (structureBuffer)
        free(structureBuffer);

	if(sequenceString.length() != structureString.length()) {
		write_log_file("Sequence and structure length do not match", fileName, 'E');
		return false;
	}

	return true;
}

bool save_file (const char *fileName, char *outputPath, const char *sequence, char *restricted, char *structure, double energy, int chosen_method) {
	FILE *ioFile;
	char filePath[MAXSLEN];
	strcpy(filePath, outputPath);
	strcat(filePath, fileName);

	std::cout << "Output File: " << filePath << '\n' << std::flush;

	ioFile = fopen(filePath, "w");
	if (ioFile != NULL) {
		fprintf(ioFile, "Sequence: %s\nInput_structure: %s\nOutput_structure: %s\nEnergy: %.2f\nMethod: %i", sequence, restricted, structure, energy, chosen_method);
	} else {
		write_log_file("Could not open file", fileName, 'E');
		return false;
	}

	// Clean up
	fclose(ioFile);
		
	return true;
} 

void write_log_file(const char *message, const char *fileName, const char option) {
	time_t now = time(0);
	struct tm tstruct = *localtime(&now);
	char buf[80];

	strftime(buf, sizeof(buf), "%c", &tstruct);
	
	switch (option) {
		case 'i':
		case 'I':
			fprintf(logFile, "[%s] Info: %s\n", buf, message);
			break;

		case 'a':
		case 'A':
			fprintf(logFile, "[%s] Action: %s in file %s\n", buf, message, fileName);
			break;

		case 'e':
		case 'E':
			fprintf(logFile, "[%s] Error: %s in file %s\n", buf, message, fileName);
			break;	
	}

	fflush(logFile);
}

/* in cases like the following
__(((____)))_____((((_____(((____)))_____(((_______)))___))))__
we have more than on B and Bp to consider
so we need to find each structure's B and Bp
the output of this function is an array, the first element of the array gives the number of substructures
then each pair gives their B and Bp 
output must be allocated first with malloc before this function is used.*/
bool find_each_substructure (char *input_structure, int begin, int end, int **output) {
	int B, Bp, ip;
	int num_structure = 0;
	int result[end-begin+1][2];
	for (auto &component : result) {
		component[0] = 0;
		component[1] = 0;
	}

	for (int i = begin; i <= end; i++) {

		if (input_structure[i] == '(') {
			if (customStack.empty()) {
				B = i;	
			}
			customStack.push(i);

		} else if (input_structure[i] == ')') {
			ip = customStack.top();
			customStack.pop();
			
			if (customStack.empty()) {
				Bp = i;
		
				if (B != ip) {
					char error[200];
					sprintf(error, "There is something wrong with finding B=%d and Bp=%d which pairs with %d\n", B, Bp, ip);
					write_log_file(error, "N/A", 'E');
					return false;
				} else {
					result[num_structure][0] = B;
					result[num_structure][1] = Bp;
					B = -1;
					Bp = -1;
					num_structure++;
				}
			}
		}

	} 

	std::copy(&result[0][0], &result[0][0] + num_structure*2, &output[0][0]);
	return true;
}

/* the structure given to this function is what is returned from simfold restricted, so has . and ( )
This function gets the beginning and end of the given structure, together with the output structure from simfold restricted when it was given the first structure as input.
it then goes back from the beginning point till it finds the start of the sequence or 3 unpaired bases, 
and moves forward from the end point till it finds the end f the sequence or 3 unpaired bases.
it returns the new beginning and end points.*/
void find_independant_structures (char *structure, int begin, int end, int *B, int *Bp) {
	*B = begin;
	*Bp = end;
	int found = 0;
	int unpaired_count = 0;
	int L_index = *B-1;
	int R_index = *Bp-1;

	while (L_index >= 0 && R_index < strlen(structure) && !found) {
		if (structure[L_index] == '(') {
			if (structure[R_index] == ')') {
				*B = L_index;
				*Bp = R_index;
				L_index--;
				R_index++;

			} else if (structure[R_index] == '.') {
				R_index++;

				while (structure[R_index] != ')' && R_index < strlen(structure) && !found) {
					R_index++;

					if ((structure[R_index] == '(') || (R_index == strlen(structure))) {
						found = 1;
					}
				
					unpaired_count++;
					if (unpaired_count == 3) {
						found = 1;
					}
				}	
			
				if (!found) {
					*B = L_index;
					*Bp = R_index;
					L_index--;
					R_index++;
				} 

			} else {
				found = 1;
			}

		} else if (structure[L_index] == '.') {
			if (structure[R_index] == ')') {
				while (structure[L_index] != '(' && L_index >= 0 && !found) {
					L_index--;

					if ((structure[L_index] == ')') || (L_index < 0)){
						found = 1;
					}
				}

				if (!found) {
					*B = L_index;
					*Bp = R_index;
					L_index--;
					R_index++;
				}

			} else if (structure[R_index] == '.') {
				R_index++;
				L_index--;
				unpaired_count++;
	
				if (unpaired_count == 3) {
					found = 1;
				}

			} else {
				found = 1;
			}
		} else {
			found = 1;
		}
	}

}

void find_sub_sequence_structure (const char *input_sequence, char *input_structure, char *output_sequence, char *output_structure, int *begin, int *end) {
	int structure_length = strlen(input_structure);	
	int sub_length = 0;
	*begin = -1;
	*end = 0;

	for (int i=0; i < structure_length; i++) {
		if ((input_structure[i] == '(') && (*begin == -1)) {
			*begin = i;
		}

		if (input_structure[i] == ')') {
			*end = i;
		}
	}

	sub_length = *end - *begin + 1;
	strncpy(output_sequence, &input_sequence[*begin], sub_length);
	strncpy(output_structure, &input_structure[*begin], sub_length);
	output_sequence[sub_length] = '\0';
	output_structure[sub_length] = '\0';
}

bool find_new_structure (char *input_structure, char *output_structure) {
	structureString.assign(input_structure, strlen(input_structure));
	structureRegex = "\\[";

	if(!std::regex_search(structureString, match, structureRegex)) {
		return false;
	}

	//substitute all ".()" to "_"
	replacementText = "_";
	structureRegex = "(\\.|\\(|\\))";
	structureString = std::regex_replace(structureString, structureRegex, replacementText);

	// substitute all [ to (
	replacementText = "(";
	structureRegex = "\\[";
	structureString = std::regex_replace(structureString, structureRegex, replacementText);

	// substitute all ] to )
	replacementText = ")";
	structureRegex = "\\]";
	structureString = std::regex_replace(structureString, structureRegex, replacementText);

	strcpy(output_structure, structureString.c_str());
	//std::cout << "New Structure String: " << structureString << '\n' << std::flush;
	return true;
}

int find_no_base_pairs (char *structure) {
	structureString.assign(structure, strlen(structure));
	return std::count(structureString.begin(), structureString.end(), '(');
}


// gets a seq and structure
// runs simfold on the sequence
// compares the number of base pairs in the simfold_seq with the input sequence
// returns the ratio of the counts;
double find_structure_sparsity (char *sequence, char *structure) {
	double ratio = 0.0;
	char simfold_sequence[MAXSLEN];
	char simfold_structure[MAXSLEN];
	double simfold_energy;
	int input_structure_basepair_count;
	int simfold_structure_basepair_count;

	call_simfold(SIMFOLD, sequence, NULL, simfold_structure, &simfold_energy);
	input_structure_basepair_count = find_no_base_pairs(structure);
	simfold_structure_basepair_count = find_no_base_pairs(simfold_structure);
	ratio = (double) input_structure_basepair_count/simfold_structure_basepair_count;

	//std::cout << "Sparsity ratio is: " << ratio << '\n' << std::flush;

	return ratio;
}

