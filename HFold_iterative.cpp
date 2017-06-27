// a simple driver for the HFold
#include <boost/regex.hpp>
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

// Hosna June 20th, 2007
//#include "W_final.h"
//#include "hfold.h"
#include "hfold_iterative.h"

#include <ctime>
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
static boost::smatch match;
static boost::regex sequenceRegex;
static boost::regex structureRegex;
static boost::regex replacedRegex;
static std::string sequenceString;
static std::string structureString;
static std::string replacementText;
static std::string energyString;
static std::stack<int> customStack;

struct thread_info {
	pthread_t thread_id;
	long      thread_num;
	char     *structure;
	char     *sequence; 
	char      method_structure[MAXSLEN];
	double    method_energy;
	int     **result;
	bool      reattempt_run;
};
static void *threadFunction(void *arg);

struct thread_info_inter {
	pthread_t thread_id;
	long      thread_num;
	char     *structure_one;
	char     *sequence_one; 
	char     *structure_two;
	char     *sequence_two;
	char      method_structure[MAXSLEN];
	double    method_energy;
	int     **result;
	bool      reattempt_run;
};
static void *threadFunctionInteracting(void *arg);

void printUsage();

void run_hfold_interacting();

void method1_calculation (char *sequence, char *structure, char *method1_structure, double *method1_energy);
void method2_calculation (char *sequence, char *structure, char *method2_structure, double *method2_energy);
void method3_calculation (char *sequence, char *structure, char *method3_structure, double *method3_energy);
void method4_calculation (char *sequence, char *structure, char *method4_structure, double *method4_energy, int **result);

void method1_calculation_interacting (char *sequence_one, char *structure_one, char *sequence_two, char *structure_two, char *method1_structure, double *method1_energy, bool reattempt_run);
void method2_calculation_interacting (char *sequence_one, char *structure_one, char *sequence_two, char *structure_two, char *method2_structure, double *method2_energy, bool reattempt_run);
void method3_calculation_interacting (char *sequence_one, char *structure_one, char *sequence_two, char *structure_two, char *method3_structure, double *method3_energy, bool reattempt_run);
void method4_calculation_interacting (char *sequence_one, char *structure_one, char *sequence_two, char *structure_two, char *method4_structure, double *method4_energy, int **result, bool reattempt_run);

bool call_HFold (char *programPath, char *input_sequence, char *input_structure, char *output_structure, double *output_energy, bool reattempt_run);
bool call_HFold_interacting (char *programPath, char *input_sequence_one, char *input_structure_one, char *input_sequence_two, char *input_structure_two, char *output_structure, double *output_energy, bool reattempt_run);
bool call_simfold (char *programPath, char *input_sequence, char *input_structure, char *output_structure, double *output_energy, bool reattempt_run);

void replace_simfold_partial_structure_with_original (char *input_structure, char *simfold_structure, char *replaced_structure, int begin, int end);
void replace_simfold_structure_with_original (char *replaced_structure, char *simfold_structure, int begin, int end);
bool get_sequence_structure (char *fileName, char *sequence, char *structure);
bool get_sequence_structure_interacting (char *fileName, char *sequence_one, char *structure_one, char *sequence_two, char *structure_two);
bool save_file (char *fileName, char *outputPath, char *sequence, char *restricted, char *structure, double energy, int chosen_method);
void write_log_file(char *message, char *file, char option);

bool find_each_substructure (char *input_structure, int begin, int end, int **output);
void find_independant_structures (char *structure, int begin, int end, int *B, int *Bp);
void find_sub_sequence_structure (char *input_sequence, char *input_structure, char *output_sequence, char *output_structure, int *begin, int *end);
bool find_new_structure (char *input_structure, char *output_structure);
int find_no_base_pairs (char *structure);
int find_files(char **files_found, char *path);
double find_structure_sparsity (char *sequence, char *structure); // Unused. keep for now.
std::string exec(const char* cmd);

void validator();
void file_combiner();
void compare_output();
void create_oligo_files();
void segfault_sigaction(int signal, siginfo_t *si, void *arg);


/*
exit(0)	success
exit(1)	general error
exit(3)	thread error
exit(4) i/o error
exit(5) pipe error
*/

/* I am changing iterative HFold so that we run 4 methods and choose the structure with the lowest energy as the winner
1) run HFold_PKonly on input structure, if pked, keep the pk bases as input structure and run HFold on them and get the result
2) run HFold on the original input and get the result
3) get the beginning and end of the given structure, give the subsequence and the substructure to simfold restricted as input, get the simfold structure, then give that as input to HFold_PKonly, get the pked bases as new input to give HFold and get the result
4) run simfold restricted with the given input sequence and structure, and get the simfold structure. Only choose part of the structure that contains the original given input structure, then give that to HFold_PKonly as input, get pked bases and run HFold on them and get the result. */
int main (int argc, char **argv) {

	//pthread_t *thread_id;
	struct thread_info *tinfo;
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
	double final_energy;

	int files_length, method_chosen;
	long threadNum;
	long numThreads = 4;
	long numRerunThreads = 0;
	 
	//int **result;
	char **files_found;
	

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

	tinfo = (thread_info *) calloc(numThreads, sizeof(struct thread_info));
	if (tinfo == NULL) {
        write_log_file("Could not calloc thread_info", "", 'I');
		exit(3);
	}

	*method1_energy = INF;
	*method2_energy = INF;
	*method3_energy = INF;
	*method4_energy = INF;
	final_energy = 0;
	numRerunThreads = 0;

	method1_structure[0] = '\0';
	method2_structure[0] = '\0';
	method3_structure[0] = '\0';
	method4_structure[0] = '\0';
	final_structure[0] = '\0';		

	//create thread
	for (threadNum = 0; threadNum < numThreads; threadNum++) {
		tinfo[threadNum].thread_num = threadNum + 1;
		tinfo[threadNum].structure = structure;
		tinfo[threadNum].sequence = sequence;
		tinfo[threadNum].method_structure[0] = '\0';
		tinfo[threadNum].method_energy = INF;
		tinfo[threadNum].reattempt_run = true;
		//tinfo[threadNum].result = result;

		if (pthread_create(&tinfo[threadNum].thread_id, NULL, &threadFunction, &tinfo[threadNum]) != 0) {
			write_log_file("Could not create thread", "", 'I');
			exit(3);
		}
	}

	//join thread
	for (threadNum = 0; threadNum < numThreads; threadNum++) {
		if (pthread_join(tinfo[threadNum].thread_id, &res) != 0){
			write_log_file("Could not join threads", "", 'I');
			exit(3);
		}
		
		//if (res)
		//	free(res);      /* Free memory allocated by thread */
	}

	//re-run thread
	for (threadNum = 0; threadNum < numThreads; threadNum++) {
		if ((strcmp(tinfo[threadNum].method_structure, "") == 0) && (tinfo[threadNum].method_energy == INF)) {
			write_log_file("Rerunning Thread", "", 'I');
			numRerunThreads++;
			tinfo[threadNum].thread_num = threadNum + 1;
			tinfo[threadNum].structure = structure;
			tinfo[threadNum].sequence = sequence;
			tinfo[threadNum].method_structure[0] = '\0';
			tinfo[threadNum].method_energy = INF;
			tinfo[threadNum].reattempt_run = false;

			if (pthread_create(&tinfo[threadNum].thread_id, NULL, &threadFunctionInteracting, &tinfo[threadNum]) != 0) {
				write_log_file("Could not create thread", "", 'I');
				exit(3);
			}

			if (pthread_join(tinfo[threadNum].thread_id, &res) != 0){
				write_log_file("Could not join threads", "", 'I');
				exit(3);
			}
		}	
	}

	/*************** Compare Methods ***************/
	final_energy = INF;
	
	for (threadNum = 0; threadNum < numThreads; threadNum++) {
		if ((strcmp(tinfo[threadNum].method_structure, "") != 0) && (tinfo[threadNum].method_energy < final_energy) && (tinfo[threadNum].method_energy != 0)) {
			final_energy = tinfo[threadNum].method_energy;
			strcpy(final_structure, tinfo[threadNum].method_structure);
			method_chosen = tinfo[threadNum].thread_num;
		}
	}

	if (final_energy == INF) {
		write_log_file("Could not find energy", "", 'E');
		final_energy = 0;
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
	
	//free(result);
	free(tinfo);
	free(method1_structure);
	free(method2_structure);
	free(method3_structure);
	free(method4_structure);
	free(method1_energy);
	free(method2_energy);
	free(method3_energy);
	free(method4_energy);

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

static void *threadFunction(void *arg) {
	struct thread_info *tinfo = (thread_info *) arg;

	switch (tinfo->thread_num) {
		case 1:
			method1_calculation(tinfo->sequence, tinfo->structure, tinfo->method_structure, &tinfo->method_energy);
			break;
		case 2:
			method2_calculation(tinfo->sequence, tinfo->structure, tinfo->method_structure, &tinfo->method_energy);
			break;
		case 3:
			method3_calculation(tinfo->sequence, tinfo->structure, tinfo->method_structure, &tinfo->method_energy);
			break;
		case 4:
			tinfo->result = (int**) malloc(sizeof(int*) * MAXSLEN);
			for (int i=0; i < MAXSLEN; i++) {
				tinfo->result[i] = (int*) malloc(sizeof(int) * 2);
			}

			method4_calculation(tinfo->sequence, tinfo->structure, tinfo->method_structure, &tinfo->method_energy, tinfo->result);

			for (int i=0; i < MAXSLEN; i++) {
				if (tinfo->result[i])
					free(tinfo->result[i]);
			}

			break;
	}
}

void method1_calculation (char *sequence, char *structure, char *method1_structure, double *method1_energy) {
	char new_input_structure[MAXSLEN] = "\0";	
	char hfold_structure[MAXSLEN] = "\0";	
	char hfold_pkonly_structure[MAXSLEN] = "\0";	
	
	double hfold_energy = 0;	
	double hfold_pkonly_energy = 0;

	bool has_pk;

	if (!call_HFold(HFOLD_PKONLY, sequence, structure, hfold_pkonly_structure, &hfold_pkonly_energy, true)) {
		*method1_energy = hfold_pkonly_energy;		
		return;
	}

	has_pk = find_new_structure(hfold_pkonly_structure, new_input_structure);
	//std::cout << "hfold_pkonly_structure = " << hfold_pkonly_structure << " new_input_structure = " << new_input_structure << " has_pk = " << has_pk << '\n' << std::flush;

	if (has_pk) {
		if (!call_HFold(HFOLD, sequence, new_input_structure, hfold_structure, &hfold_energy, true)) {
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
	
	if (!call_HFold(HFOLD, sequence, structure, hfold_structure, &hfold_energy, true)) {
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
	//std::cout << "found sub sequence = " << sub_sequence << " and sub structure = " << sub_structure << "  with begin = " << begin << " and end = " << end << '\n' << std::flush;

	if (!call_simfold(SIMFOLD, sub_sequence, sub_structure, simfold_structure, &simfold_energy, true)) {
		*method3_energy = simfold_energy;		
		return;
	}
	replace_simfold_partial_structure_with_original(structure, simfold_structure, replaced_structure, begin, end);
	//std::cout << "the replaced string is: " << replaced_structure << '\n' << std::flush;

	if (!call_HFold(HFOLD_PKONLY, sequence, replaced_structure, hfold_pkonly_structure, &hfold_pkonly_energy, true)) {
		*method3_energy = hfold_pkonly_energy;
		return;
	}
	has_pk = find_new_structure(hfold_pkonly_structure, new_input_structure);
	//std::cout << "hfold_pkonly_structure = " << hfold_pkonly_structure << " new_input_structure = " << new_input_structure << " has_pk = " << has_pk << '\n' << std::flush;

	if (has_pk) {
		if (!call_HFold(HFOLD, sequence, new_input_structure, hfold_structure, &hfold_energy, true)) {
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
	
	if (!call_simfold(SIMFOLD, sequence, structure, simfold_structure, &simfold_energy, true)) {
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

	if (!call_HFold(HFOLD_PKONLY, sequence, replaced_structure, hfold_pkonly_structure, &hfold_pkonly_energy, true)) {
		*method4_energy = hfold_pkonly_energy;
		return;
	}
	has_pk = find_new_structure(hfold_pkonly_structure, new_input_structure);
	//std::cout << "hfold_pkonly_structure = " << hfold_pkonly_structure << " new_input_structure = " << new_input_structure << " has_pk = " << has_pk << '\n' << std::flush;

	if (has_pk) {
		if (!call_HFold(HFOLD, sequence, new_input_structure, hfold_structure, &hfold_energy, true)) {
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

void run_hfold_interacting () {
	struct thread_info_inter *tinfo;
	void *res;
	
	char *output_path;
	char sequence_combined[MAXSLEN];
	char sequence_one[MAXSLEN];
	char sequence_two[MAXSLEN];
	char structure_combined[MAXSLEN];
	char structure_one[MAXSLEN];
	char structure_two[MAXSLEN];
	char method1_structure[MAXSLEN];
	char method2_structure[MAXSLEN];
	char method3_structure[MAXSLEN];
	char method4_structure[MAXSLEN];
	char final_structure[MAXSLEN];

	double method1_energy;
	double method2_energy;
	double method3_energy;
	double method4_energy;
	double final_energy;

	int result_length, files_length, method_chosen;
	long threadNum;
	long numThreads = 4;
	long numRerunThreads = 0;

	//int **result;
	char **files_found;

	logFile = fopen(LOGFILEPATH, "w");
	if (logFile == NULL) {
		giveup("Could not open logfile", LOGFILEPATH);
	}
	
	output_path = (char*) malloc(sizeof(char) * 1000);
	file = (char*) malloc(sizeof(char) * 1000);
	files_found = (char**) malloc(MAXSLEN*sizeof(char*));
	//result = (int**) malloc(MAXSLEN*sizeof(int*));
	for (int i=0; i < MAXSLEN; i++) {
		files_found[i] = (char*) malloc(sizeof(char) * 255);
		//result[i] = (int*) malloc(sizeof(int) * 2);
	}
	
	write_log_file("Starting Program", "", 'I');

	tinfo = (thread_info_inter *) calloc(numThreads, sizeof(struct thread_info_inter));
	if (tinfo == NULL) {
        write_log_file("Could not calloc thread_info", "", 'I');
		exit(3);
	}

	files_length = find_files(files_found, INPUTPATH);

	for (int k = 0; k < files_length; k++) {
		strcpy(output_path, OUTPUTPATH);
		strcat(output_path, files_found[k]);

		if(access(output_path, F_OK) != -1) {
			continue;
		}	
			
		strcpy(file, INPUTPATH);
		strcat(file, files_found[k]);

		numRerunThreads = 0;
		method1_energy = 0;
		method2_energy = 0;
		method3_energy = 0;
		method4_energy = 0;
		final_energy = 0;

		sequence_combined[0] = '\0';
		sequence_one[0] = '\0';
		sequence_two[0] = '\0';	
		structure_combined[0] = '\0';
		structure_one[0] = '\0';	
		structure_two[0] = '\0';	
		method1_structure[0] = '\0';
		method2_structure[0] = '\0';
		method3_structure[0] = '\0';
		method4_structure[0] = '\0';
		final_structure[0] = '\0';

		if (!get_sequence_structure_interacting (file, sequence_one, structure_one, sequence_two, structure_two)) {
			continue;
		}		

		for (threadNum = 0; threadNum < numThreads; threadNum++) {
			tinfo[threadNum].thread_num = threadNum + 1;
			tinfo[threadNum].structure_one = structure_one;
			tinfo[threadNum].sequence_one = sequence_one;
			tinfo[threadNum].structure_two = structure_two;
			tinfo[threadNum].sequence_two = sequence_two;
			tinfo[threadNum].method_structure[0] = '\0';
			tinfo[threadNum].method_energy = INF;
			tinfo[threadNum].reattempt_run = true;
			//tinfo[threadNum].result = result;

		    if (pthread_create(&tinfo[threadNum].thread_id, NULL, &threadFunctionInteracting, &tinfo[threadNum]) != 0) {
		        write_log_file("Could not create thread", "", 'I');
				exit(3);
			}
		}

		for (threadNum = 0; threadNum < numThreads; threadNum++) {
			if (pthread_join(tinfo[threadNum].thread_id, &res) != 0){
		        write_log_file("Could not join threads", "", 'I');
				exit(3);
			}
			
			//if (res)
			//	free(res);      /* Free memory allocated by thread */
		}

		for (threadNum = 0; threadNum < numThreads; threadNum++) {
			if ((strcmp(tinfo[threadNum].method_structure, "") == 0) && (tinfo[threadNum].method_energy > INF)) {
				write_log_file("Reruning thread", "", 'I');
				numRerunThreads++;
				tinfo[threadNum].thread_num = threadNum + 1;
				tinfo[threadNum].structure_one = structure_one;
				tinfo[threadNum].sequence_one = sequence_one;
				tinfo[threadNum].structure_two = structure_two;
				tinfo[threadNum].sequence_two = sequence_two;
				tinfo[threadNum].method_structure[0] = '\0';
				tinfo[threadNum].method_energy = INF;
				tinfo[threadNum].reattempt_run = false;

				if (pthread_create(&tinfo[threadNum].thread_id, NULL, &threadFunctionInteracting, &tinfo[threadNum]) != 0) {
				    write_log_file("Could not create thread", "", 'I');
					exit(3);
				}

				if (pthread_join(tinfo[threadNum].thread_id, &res) != 0){
				    write_log_file("Could not join threads", "", 'I');
					exit(3);
				}
			}	
        }

		/*************** First Method ***************/
		//method1_calculation_interacting(sequence_one, structure_one, sequence_two, structure_two, method1_structure, &method1_energy);

		/*************** Second Method ***************/
		//method2_calculation_interacting(sequence_one, structure_one, sequence_two, structure_two, method2_structure, &method2_energy);
	
		/*************** Third Method ***************/
		//method3_calculation_interacting(sequence_one, structure_one, sequence_two, structure_two, method3_structure, &method3_energy);

		/*************** Fourth Method ***************/
		//method4_calculation_interacting(sequence_one, structure_one, sequence_two, structure_two, method4_structure, &method4_energy, result);

		/*************** Compare Methods ***************/
		final_energy = INF;
	
		for (threadNum = 0; threadNum < numThreads; threadNum++) {
			if ((strcmp(tinfo[threadNum].method_structure, "") != 0) && (tinfo[threadNum].method_energy < final_energy) && (tinfo[threadNum].method_energy != 0)) {
				final_energy = tinfo[threadNum].method_energy;
				strcpy(final_structure, tinfo[threadNum].method_structure);
				method_chosen = tinfo[threadNum].thread_num;
			}
		}		

		/*if (method1_energy < final_energy && method1_energy != 0) {
			final_energy = method1_energy;
			strcpy(final_structure, method1_structure);
			method_chosen = 1;
		}

		if (method2_energy < final_energy && method2_energy != 0) {
			final_energy = method2_energy;	
			strcpy(final_structure, method2_structure);		
			method_chosen = 2;
		}

		if (method3_energy < final_energy && method3_energy != 0) {
			final_energy = method3_energy;	
			strcpy(final_structure, method3_structure);		
			method_chosen = 3;
		}

		if (method4_energy < final_energy && method4_energy != 0) {
			final_energy = method4_energy;	
			strcpy(final_structure, method4_structure);		
			method_chosen = 4;
		}*/

		if (final_energy == INF) {
			write_log_file("Could not find energy", "", 'E');
			final_energy = 0;
		}

		strcpy(sequence_combined, sequence_one);
		strcat(sequence_combined, linker);
		strcat(sequence_combined, sequence_two);

		strcpy(structure_combined, structure_one);
		strcat(structure_combined, ".....");
		strcat(structure_combined, structure_two);

		std::cout << "final structure: " << final_structure << " final energy: " << final_energy << " method chosen: " << method_chosen << '\n' << std::flush;
		save_file(files_found[k], OUTPUTPATH, sequence_combined, structure_combined, final_structure, final_energy, method_chosen);
	}

	// Clean up
	

	for (int i=0; i < MAXSLEN; i++) {
		if (files_found[i])
			free(files_found[i]);
		//if (result[i])
			//free(result[i]);
	}
	free(files_found);
	//free(result);
	free(file);
	free(tinfo);

	// This will cause a seg fault if you have the log file open and are watching it.
	if (logFile) {
		fclose(logFile);
	}
}

static void *threadFunctionInteracting(void *arg) {
	struct thread_info_inter *tinfo = (thread_info_inter *) arg;

	switch (tinfo->thread_num) {
		case 1:
			method1_calculation_interacting(tinfo->sequence_one, tinfo->structure_one, tinfo->sequence_two, tinfo->structure_two, tinfo->method_structure, &tinfo->method_energy, tinfo->reattempt_run);
			break;
		case 2:
			method2_calculation_interacting(tinfo->sequence_one, tinfo->structure_one, tinfo->sequence_two, tinfo->structure_two, tinfo->method_structure, &tinfo->method_energy, tinfo->reattempt_run);
			break;
		case 3:
			method3_calculation_interacting(tinfo->sequence_one, tinfo->structure_one, tinfo->sequence_two, tinfo->structure_two, tinfo->method_structure, &tinfo->method_energy, tinfo->reattempt_run);
			break;
		case 4:
			tinfo->result = (int**) malloc(sizeof(int*) * MAXSLEN);
			for (int i=0; i < MAXSLEN; i++) {
				tinfo->result[i] = (int*) malloc(sizeof(int) * 2);
			}

			method4_calculation_interacting(tinfo->sequence_one, tinfo->structure_one, tinfo->sequence_two, tinfo->structure_two, tinfo->method_structure, &tinfo->method_energy, tinfo->result, tinfo->reattempt_run);
			
			for (int i=0; i < MAXSLEN; i++) {
				if (tinfo->result[i])
					free(tinfo->result[i]);
			}	

			break;
	}
}

void method1_calculation_interacting (char *sequence_one, char *structure_one, char *sequence_two, char *structure_two, char *method1_structure, double *method1_energy, bool reattempt_run) {
	char new_input_structure[MAXSLEN] = "\0";	
	char new_input_structure_one[MAXSLEN] = "\0";
	char new_input_structure_two[MAXSLEN] = "\0";	
	char hfold_structure[MAXSLEN] = "\0";	
	char hfold_pkonly_structure[MAXSLEN] = "\0";	
	
	double hfold_energy = 0;	
	double hfold_pkonly_energy = 0;
	bool has_pk;

	if (!call_HFold_interacting(HFOLD_INTERACTING_PKONLY, sequence_one, structure_one, sequence_two, structure_two, hfold_pkonly_structure, &hfold_pkonly_energy, reattempt_run)) {
		*method1_energy = hfold_pkonly_energy;
		return;
	}
	has_pk = find_new_structure(hfold_pkonly_structure, new_input_structure);

	strncpy(new_input_structure_one, new_input_structure, strlen(structure_one));
	strncpy(new_input_structure_two, &new_input_structure[strlen(structure_one)+linker_length], strlen(structure_two));

	//std::cout << "hfold_pkonly_structure = " << hfold_pkonly_structure << " new_input_structure = " << new_input_structure << " has_pk = " << has_pk << '\n' << std::flush;

	if (has_pk) {
		if (!call_HFold_interacting(HFOLD_INTERACTING, sequence_one, new_input_structure_one, sequence_two, new_input_structure_two, hfold_structure, &hfold_energy, reattempt_run)) {
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

void method2_calculation_interacting (char *sequence_one, char *structure_one, char *sequence_two, char *structure_two, char *method2_structure, double *method2_energy, bool reattempt_run) {
	char hfold_structure[MAXSLEN] = "\0";	
	
	double hfold_energy = 0;	
	
	if (!call_HFold_interacting(HFOLD_INTERACTING, sequence_one, structure_one, sequence_two, structure_two, hfold_structure, &hfold_energy, reattempt_run)) {
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

void method3_calculation_interacting (char *sequence_one, char *structure_one, char *sequence_two, char *structure_two, char *method3_structure, double *method3_energy, bool reattempt_run) {
	char sub_sequence[MAXSLEN] = "\0";
	char sub_structure[MAXSLEN] = "\0";
	char sequence_combined[MAXSLEN] = "\0";
	char structure_combined[MAXSLEN] = "\0";

	char replaced_structure[MAXSLEN] = "\0";
	char replaced_structure_one[MAXSLEN] = "\0";
	char replaced_structure_two[MAXSLEN] = "\0";

	char new_input_structure[MAXSLEN] = "\0";
	char new_input_structure_one[MAXSLEN] = "\0";
	char new_input_structure_two[MAXSLEN] = "\0";

	char hfold_structure[MAXSLEN] = "\0";
	char hfold_pkonly_structure[MAXSLEN] = "\0";
	char simfold_structure[MAXSLEN] = "\0";

	double hfold_energy = 0;	
	double hfold_pkonly_energy = 0;
	double simfold_energy = 0;
	
	int begin = -1, end = 0;
	bool has_pk;

	strcpy(sequence_combined, sequence_one);
	strcat(sequence_combined, linker);
	strcat(sequence_combined, sequence_two);

	strcpy(structure_combined, structure_one);
	strcat(structure_combined, ".....");
	strcat(structure_combined, structure_two);

	find_sub_sequence_structure(sequence_combined, structure_combined, sub_sequence, sub_structure, &begin, &end);
	//std::cout << "found sub sequence = " << sub_sequence << " and sub structure = " << sub_structure << "  with begin = " << begin << " and end = " << end << '\n' << std::flush;

	if (!call_simfold(SIMFOLD, sub_sequence, sub_structure, simfold_structure, &simfold_energy, reattempt_run)) {
		*method3_energy = simfold_energy;
		return;
	}
	replace_simfold_partial_structure_with_original(structure_combined, simfold_structure, replaced_structure, begin, end);

	strncpy(replaced_structure_one, replaced_structure, strlen(structure_one));
	strncpy(replaced_structure_two, &replaced_structure[strlen(structure_one)+linker_length], strlen(structure_two));

	//std::cout << "the replaced string is: " << replaced_structure << '\n' << std::flush;
	if (!call_HFold_interacting(HFOLD_INTERACTING_PKONLY, sequence_one, replaced_structure_one, sequence_two, replaced_structure_two, hfold_pkonly_structure, &hfold_pkonly_energy, reattempt_run)) {
		*method3_energy = hfold_pkonly_energy;
		return;
	}

	has_pk = find_new_structure(hfold_pkonly_structure, new_input_structure);

	strncpy(new_input_structure_one, new_input_structure, strlen(structure_one));
	strncpy(new_input_structure_two, &new_input_structure[strlen(structure_one)+linker_length], strlen(structure_two));

	//std::cout << "hfold_pkonly_structure = " << hfold_pkonly_structure << " new_input_structure = " << new_input_structure << " has_pk = " << has_pk << '\n' << std::flush;

	if (has_pk) {
		if (!call_HFold_interacting(HFOLD_INTERACTING, sequence_one, new_input_structure_one, sequence_two, new_input_structure_two, hfold_structure, &hfold_energy, reattempt_run)) {
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

void method4_calculation_interacting (char *sequence_one, char *structure_one, char *sequence_two, char *structure_two, char *method4_structure, double *method4_energy, int **result, bool reattempt_run) {
	// we know the beginning and end of the given input substructure in the structure
	// now we need to run simfold restricted with the given structure, find the stem which contains the given substructure
	// then only include this structure as input structure, and do as we did in method 3
	char sub_sequence[MAXSLEN] = "\0";
	char sub_structure[MAXSLEN] = "\0";
	char sequence_combined[MAXSLEN] = "\0";
	char structure_combined[MAXSLEN] = "\0";

	char replaced_structure[MAXSLEN] = "\0";
	char replaced_structure_one[MAXSLEN] = "\0";
	char replaced_structure_two[MAXSLEN] = "\0";

	char new_input_structure[MAXSLEN] = "\0";
	char new_input_structure_one[MAXSLEN] = "\0";
	char new_input_structure_two[MAXSLEN] = "\0";

	char hfold_structure[MAXSLEN] = "\0";
	char hfold_pkonly_structure[MAXSLEN] = "\0";
	char simfold_structure[MAXSLEN] = "\0";

	double hfold_energy = 0;	
	double hfold_pkonly_energy = 0;
	double simfold_energy = 0;
	
	int begin = -1, end = 0, B = 0, Bp = 0;
	int result_length = 0;
	bool has_pk;	

	strcpy(sequence_combined, sequence_one);
	strcat(sequence_combined, linker);
	strcat(sequence_combined, sequence_two);
	
	strcpy(structure_combined, structure_one);
	strcat(structure_combined, ".....");
	strcat(structure_combined, structure_two);

	strcpy(structure_combined, structure_combined);

	find_sub_sequence_structure(sequence_combined, structure_combined, sub_sequence, sub_structure, &begin, &end);

	if (!call_simfold(SIMFOLD, sequence_combined, structure_combined, simfold_structure, &simfold_energy, reattempt_run)) {
		*method4_energy = simfold_energy;
		return;
	}
	
	if (!find_each_substructure(structure_combined, begin, end, result)) {
		return;
	}

	result_length = sizeof(result) / sizeof(result[0]);
	//strcpy(replaced_structure, structure_combined);
	for (int i = 0; i < result_length; i++) {
		find_independant_structures (simfold_structure, result[i][0], result[i][1], &B, &Bp);
		//std::cout << "found B = " << B << " and Bp = " << Bp << '\n' << std::flush;
		replace_simfold_structure_with_original(replaced_structure, simfold_structure, B, Bp);
	}

	strncpy(replaced_structure_one, replaced_structure, strlen(structure_one));
	strncpy(replaced_structure_two, &replaced_structure[strlen(structure_one)+linker_length], strlen(structure_two));

	//std::cout << "the replaced structure is: " << replaced_structure << '\n' << std::flush;
	if (!call_HFold_interacting(HFOLD_INTERACTING_PKONLY, sequence_one, replaced_structure_one, sequence_two, replaced_structure_two, hfold_pkonly_structure, &hfold_pkonly_energy, reattempt_run)) {
		*method4_energy = hfold_pkonly_energy;
		return;
	}
	has_pk = find_new_structure(hfold_pkonly_structure, new_input_structure);

	strncpy(new_input_structure_one, new_input_structure, strlen(structure_one));
	strncpy(new_input_structure_two, &new_input_structure[strlen(structure_one)+linker_length], strlen(structure_two));

	//std::cout << "hfold_pkonly_structure = " << hfold_pkonly_structure << " new_input_structure = " << new_input_structure << " has_pk = " << has_pk << '\n' << std::flush;

	if (has_pk) {
		if (!call_HFold_interacting(HFOLD_INTERACTING, sequence_one, new_input_structure_one, sequence_two, new_input_structure_two, hfold_structure, &hfold_energy, reattempt_run)) {
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

bool call_HFold (char *programPath, char *input_sequence, char *input_structure, char *output_structure, double *output_energy, bool reattempt_run) {
	int status = 0;
	int pipefd[2];
	pipe(pipefd);
	std::string result = "";
	//int stdout_bk = dup(fileno(stdout));

	pid_t pid = fork(); /* Create a child process */
	
	switch (pid) {
		case -1: /* Error */
			write_log_file("Fork failed", file, 'E');
			return false;		    

		case 0: /* Child process */
			char functionCall[10000];
			char result_array[10000];
			close(pipefd[0]);	

			//std::cout << "Calling: " << programPath << " " << input_sequence << " " << input_structure << '\n' << std::flush;
	
			// Create a pipe to get the result through stdout
			//dup2(pipefd[1], fileno(stdout));
			//execl(programPath, "./", input_sequence, input_structure, NULL);
			
			// Run HFold_pkonly
			sprintf(functionCall, "%s %s '%s'", programPath, input_sequence, input_structure);
			result += exec(functionCall);
			strcpy(result_array, result.c_str());

			if (!result.empty()){
				write(pipefd[1], result_array, 10000); 
				exit(5);
			} else {
				write_log_file("popen() failed!", file, 'E');
				close(pipefd[1]);
				exit(5);
			}

			//char error[200];
			//sprintf(error, "System call failed with error: %s", strerror(errno));
			//write_log_file(error, programPath, 'E');
		    //exit(2);

		default: /* Parent process */
			char buf[10001];
			char *token[4];
            char sequenceTest[2000];   
            char restrictedTest[2000];
			char *saveptr1, *saveptr2;
			pid_t wpid;
			int waittime = 0;
			//structureRegex = ("RES:\\s+([^\\s+]*)\\s+([^\\n+]*)");

			// Make child process the leader of its own process group. This allows
   			// signals to also be delivered to processes forked by the child process.
   			setpgid(pid, 0);

		    //waitpid(pid, &status, 0); /* Wait for the process to complete */
			//sprintf(error, "The pid of the child is: %d", pid);
			//write_log_file(error, "", 'I');
			do {
				wpid = waitpid(pid, &status, WNOHANG);
				if (wpid == 0) {
				    if (waittime < timeout) {
						//sprintf(error, "Parent waiting %d second(s).", waittime);
				        //write_log_file(error, "", 'I');
				        sleep(1);
				        waittime ++;
				    }
				    else {
						write_log_file("Killing child process", file, 'E');
				        //kill(pid, SIGKILL); 

						kill((-1*pid), SIGTERM); 	
						bool died = false;
						for (int loop = 0; !died && loop < 5; loop++) {
							sleep(1);
							if (waitpid(pid, &status, WNOHANG) == pid) died = true;
						}

						if (!died) {
							kill((-1*pid), SIGKILL);		
						}					

						//dup2(stdout_bk, fileno(stdout));
						close(pipefd[1]);
						close(pipefd[0]);
						//close(stdout_bk);
						
						if (reattempt_run) {
							write_log_file("Reattempting to run child", "", 'I');
							*output_energy = INF + 1;
							return false;
						} else {
							write_log_file("Child will not run", "", 'I');
							return false;
						}
				    }
				}
			} while (wpid == 0 && waittime <= timeout);

			if (WEXITSTATUS(status) == 2) { 
				char error[200];
				sprintf(error, "status = %d", WEXITSTATUS(status));
				write_log_file(error, file, 'E');

				//dup2(stdout_bk, fileno(stdout));
				close(pipefd[1]);
				close(pipefd[0]);
				//close(stdout_bk);
				return false;
			}

			//restore
			//fflush(stdout);
			close(pipefd[1]);
			//dup2(stdout_bk, fileno(stdout));

			// Read result from HFold_pkonly
   			read(pipefd[0], buf, 10000);
			close(pipefd[0]);
			//close(stdout_bk);

			token[0] = strtok_r(buf, "\n", &saveptr1); //seq
            token[1] = strtok_r(saveptr1, "\n", &saveptr2); //struct

            if(token[0] != '\0')
				strcpy(sequenceTest, token[0]);
            if(token[1] != '\0')
				strcpy(restrictedTest, token[1]);

            token[2] = strtok_r(sequenceTest, "  ", &saveptr1);
            token[3] = strtok_r(restrictedTest, "  ", &saveptr2);
           
            // Error messages can be in token [0] where the program used printf once then exited.
            // A check is added here to ensure that this error is handled correctly.
            //if (token[0][0] == 'S' && token[1] != NULL) {
           
            if (token[1] != '\0' && (strcmp(token[2], "Seq:") == 0) && (strcmp(token[3], "RES:") == 0)) {
                //structureString.assign(token[1], strlen(token[1]));
                token[0] = strtok_r(token[1], "  ", &saveptr1);
                token[2] = strtok_r(saveptr1, "  ", &saveptr2);
                strcpy(output_structure, token[2]);       

                token[2] = strtok_r(saveptr2, "\n", &saveptr2);
                energyString.assign(token[2], strlen(token[2]));
                *output_energy = stod(energyString);
            } else {
                write_log_file(token[0], file, 'E');
                write_log_file(token[1], file, 'E');
                return false;
            }

			// Find the structure
			/*if(boost::regex_search(structureString, match, structureRegex)) {
				structureString = match[1];
			} else {
				write_log_file(token[0], programPath, 'E');
				return false;
			}
			
			token[0] = strtok_r(token[1], "  ", &saveptr1);
			token[2] = strtok_r(saveptr1, "  ", &saveptr2); 
			token[2] = strtok_r(saveptr2, "\n", &saveptr2); 
			
			energyString.assign(token[2], strlen(token[2]));*/

			//Print result
			//std::cout << "Result: " << structureString << " " << token[2] << '\n' << std::flush;
	
			// Copy strings
			//strcpy(output_structure, structureString.c_str());
			//*output_energy = stod(energyString);
			//write_log_file("hfold struct", output_structure, 'E');
			return true;
	}
}

bool call_HFold_interacting (char *programPath, char *input_sequence_one, char *input_structure_one, char *input_sequence_two, char *input_structure_two, char *output_structure, double *output_energy, bool reattempt_run) {
	int status = 0;	
	int pipefd[2];
	pipe(pipefd);
	std::string result = "";
	//int stdout_bk = dup(fileno(stdout));

	pid_t pid = fork(); /* Create a child process */
	
	switch (pid) {
		case -1: /* Error */
		    write_log_file("Fork failed", file, 'E');
			return false;		    

		case 0: /* Child process */
			char functionCall[10000];
			char result_array[10000];
			close(pipefd[0]);

			//std::cout << "Calling: " << programPath << " " << input_sequence_one << " " << input_structure_one << " " << input_sequence_two << " " << input_structure_two << '\n' << std::flush;
	
			// Create a pipe to get the result through stdout
			//dup2(pipefd[1], fileno(stdout));
			
			//execl(programPath, "./", input_sequence_one, input_structure_one, input_sequence_two, input_structure_two, NULL);
	
			// Run HFold_pkonly
			sprintf(functionCall, "%s %s '%s' %s '%s'", programPath, input_sequence_one, input_structure_one, input_sequence_two, input_structure_two);
			result += exec(functionCall);
			strcpy(result_array, result.c_str());

			if (!result.empty()){
				write(pipefd[1], result_array, 10000); 
				exit(5);
			} else {
				write_log_file("popen() failed!", file, 'E');
				close(pipefd[1]);
				exit(5);
			}
			
			//char error[200];
			//sprintf(error, "System call failed with error: %s", strerror(errno));
			//write_log_file(error, file, 'E');
		    //exit(2);

		default: /* Parent process */
			char buf[10001];
			char *token[4];
            char sequenceTest[2000];   
            char restrictedTest[2000];
			char *saveptr1, *saveptr2;
			pid_t wpid;
			int waittime = 0;

			// Make child process the leader of its own process group. This allows
   			// signals to also be delivered to processes forked by the child process.
   			setpgid(pid, 0);

			//structureRegex = "RES:\\s+([^\\s+]*)\\s+(-?\\d{1,4}(\\.\\d{1,2}))";
		    //waitpid(pid, &status, 0); /* Wait for the process to complete */
			//sprintf(error, "The pid of the child is: %d", pid);
			//write_log_file(error, "", 'I');
			do {
				wpid = waitpid(pid, &status, WNOHANG);
				if (wpid == 0) {
				    if (waittime < timeout) {
						//sprintf(error, "Parent waiting %d second(s).", waittime);
				        //write_log_file(error, "", 'I');
				        sleep(1);
				        waittime ++;
				    }
				    else {
						write_log_file("Killing child process", file, 'E');
				        //kill(pid, SIGKILL); 

						kill((-1*pid), SIGTERM); 	
						bool died = false;
						for (int loop = 0; !died && loop < 5; loop++) {
							sleep(1);
							if (waitpid(pid, &status, WNOHANG) == pid) died = true;
						}

						if (!died) {
							kill((-1*pid), SIGKILL);		
						}					

						//dup2(stdout_bk, fileno(stdout));
						close(pipefd[1]);
						close(pipefd[0]);
						//close(stdout_bk);

						if (reattempt_run) {
							write_log_file("Reattempting to run child", "", 'I');
							*output_energy = INF + 1;
							return false;
						} else {
							write_log_file("Child will not run", "", 'I');
							return false;
						}
				    }
				}
			} while (wpid == 0 && waittime <= timeout);

			if (WEXITSTATUS(status) == 2) { 
				char error[200];
				sprintf(error, "status = %d", WEXITSTATUS(status));
				write_log_file(error, file, 'E');

				//dup2(stdout_bk, fileno(stdout));
				close(pipefd[1]);
				close(pipefd[0]);
				//close(stdout_bk);
				return false;
			}

			//restore
			//fflush(stdout);
			close(pipefd[1]);
			//dup2(stdout_bk, fileno(stdout));

			// Read result from HFold_pkonly
   			read(pipefd[0], buf, 10000);
			close(pipefd[0]);
			//close(stdout_bk);

			//strcpy(buf, result.c_str());

			token[0] = strtok_r(buf, "\n", &saveptr1); //seq
            token[1] = strtok_r(saveptr1, "\n", &saveptr2); //struct

            if(token[0] != '\0')
				strcpy(sequenceTest, token[0]);
            if(token[1] != '\0')
				strcpy(restrictedTest, token[1]);

            token[2] = strtok_r(sequenceTest, "  ", &saveptr1);
            token[3] = strtok_r(restrictedTest, "  ", &saveptr2);
           
            // Error messages can be in token [0] where the program used printf once then exited.
            // A check is added here to ensure that this error is handled correctly.
            //if (token[0][0] == 'S' && token[1] != NULL) {
           
            if (token[1] != '\0' && (strcmp(token[2], "Seq:") == 0) && (strcmp(token[3], "RES:") == 0)) {
                //structureString.assign(token[1], strlen(token[1]));
                token[0] = strtok_r(token[1], "  ", &saveptr1);
                token[2] = strtok_r(saveptr1, "  ", &saveptr2);
                strcpy(output_structure, token[2]);       

                token[2] = strtok_r(saveptr2, "\n", &saveptr2);
                energyString.assign(token[2], strlen(token[2]));
                *output_energy = stod(energyString);
            } else {
                write_log_file(token[0], file, 'E');
                write_log_file(token[1], file, 'E');
                return false;
            }

			// Error messages can be in token [0] where the program used printf once then exited.
			// A check is added here to ensure that this error is handled correctly.
			/*if (token[1] != NULL) {
				structureString.assign(token[1], strlen(token[1]));
			} else {
				write_log_file("Token[1] is null", file, 'I'); // probably remove
				write_log_file(token[0], file, 'E');
				return false;
			}

			// Find the structure
			if(boost::regex_search(structureString, match, structureRegex)) {
				structureString = match[1];
			} else {
				write_log_file(token[0], programPath, 'E');
				return false;
			}
			
			token[0] = strtok(token[1], "  ");
			token[2] = strtok(NULL, "  "); 
			token[2] = strtok(NULL, "\n"); 
			
			energyString.assign(token[2], strlen(token[2]));

			//Print result
			std::cout << "Result: " << structureString << " " << energyString << '\n' << std::flush;
	
			// Copy strings
			strcpy(output_structure, structureString.c_str());
			*output_energy = stod(energyString);*/
			return true;
	}
}

bool call_simfold (char *programPath, char *input_sequence, char *input_structure, char *output_structure, double *output_energy, bool reattempt_run) {
	int status = 0;
	int pipefd[2];
	pipe(pipefd);
	std::string result = "";
	//int stdout_bk = dup(fileno(stdout));

	pid_t pid = fork(); /* Create a child process */
	
	switch (pid) {
		case -1: /* Error */
		    write_log_file("Fork failed", file, 'E');
			return false;		    

		case 0: /* Child process */
			char functionCall[10000];
			char result_array[10000];
			close(pipefd[0]);

			// Child changes directory so the parent doesnt have to move back.
			if (chdir(programPath) == -1) {
				char error[200];
				sprintf(error, "Could not find directory with error: %s", strerror(errno));
				write_log_file(error, file, 'E');
				exit(4);
			}

			// Create function call
			strcpy(functionCall, SIMFOLD);
			strcat(functionCall, " -s ");
			strcat(functionCall, input_sequence);
	
			if (input_structure != NULL) {
				strcat(functionCall, " -r \"");
				strcat(functionCall, input_structure);
				strcat(functionCall, "\"");
			}

			result += exec(functionCall);
			strcpy(result_array, result.c_str());

			if (!result.empty()){
				write(pipefd[1], result_array, 10000); 
				exit(5);
			} else {
				write_log_file("popen() failed!", file, 'E');
				close(pipefd[1]);
				exit(5);
			}

			//std::cout << "Calling: " << functionCall << '\n' << std::flush;

			// Create a pipe to get the result through stdout
			//dup2(pipefd[1], fileno(stdout));
			//execl(SIMFOLD, "./", "-s", input_sequence, "-r", input_structure, NULL);

			//char error[200];
			//sprintf(error, "System call failed with error: %s", strerror(errno));
			//write_log_file(error, programPath, 'E');
		    //exit(1);

		default: /* Parent process */
			char buf[10001];
			char *token[4];
            char sequenceTest[2000];   
            char restrictedTest[2000];
			char *saveptr1, *saveptr2;
			pid_t wpid;
			int waittime = 0;

			/*if (input_structure == NULL) {
				structureRegex = "MFE:\\s+([^\\s+]*)\\s+([^\\n+]*)";
			} else {
				structureRegex = "RES:\\s+([^\\s+]*)\\s+([^\\n+]*)";
			}*/ 

			// Make child process the leader of its own process group. This allows
   			// signals to also be delivered to processes forked by the child process.
   			setpgid(pid, 0);

		    //waitpid(pid, &status, 0); /* Wait for the process to complete */
			//sprintf(error, "The pid of the child is: %d", pid);
			//write_log_file(error, "", 'I');
			do {
				wpid = waitpid(pid, &status, WNOHANG);
				if (wpid == 0) {
				    if (waittime < timeout) {
						//sprintf(error, "Parent waiting %d second(s).", waittime);
				        //write_log_file(error, "", 'I');
				        sleep(1);
				        waittime ++;
				    }
				    else {
						write_log_file("Killing child process", file, 'E');
				        //kill(pid, SIGKILL); 

						kill((-1*pid), SIGTERM); 	
						bool died = false;
						for (int loop = 0; !died && loop < 5; loop++) {
							sleep(1);
							if (waitpid(pid, &status, WNOHANG) == pid) died = true;
						}

						if (!died) {
							kill((-1*pid), SIGKILL);		
						}					

						//dup2(stdout_bk, fileno(stdout));
						close(pipefd[1]);
						close(pipefd[0]);
						//close(stdout_bk);
						
						if (reattempt_run) {
							write_log_file("Reattempting to run child", "", 'I');
							*output_energy = INF + 1;
							return false;
						} else {
							write_log_file("Child will not run", "", 'I');
							return false;
						}
				    }
				}
			} while (wpid == 0 && waittime <= timeout);

			if (WEXITSTATUS(status) == 2) { 
				char error[200];
				sprintf(error, "status = %d", WEXITSTATUS(status));
				write_log_file(error, file, 'E');

				//dup2(stdout_bk, fileno(stdout));
				close(pipefd[1]);
				close(pipefd[0]);
				//close(stdout_bk);
				return false;
			}
		
			//restore
			//fflush(stdout);
			close(pipefd[1]);
			//dup2(stdout_bk, fileno(stdout));

			// Read result from HFold_pkonly
   			read(pipefd[0], buf, 10000);
			close(pipefd[0]);
			//close(stdout_bk);

			token[0] = strtok_r(buf, "\n", &saveptr1); //seq
            token[1] = strtok_r(saveptr1, "\n", &saveptr2); //struct

            if(token[0] != '\0')
				strcpy(sequenceTest, token[0]);
            if(token[1] != '\0')
				strcpy(restrictedTest, token[1]);

            token[2] = strtok_r(sequenceTest, "  ", &saveptr1);
            token[3] = strtok_r(restrictedTest, "  ", &saveptr2);
           
            // Error messages can be in token [0] where the program used printf once then exited.
            // A check is added here to ensure that this error is handled correctly.
            //if (token[0][0] == 'S' && token[1] != NULL) {
           
            if (token[1] != '\0' && (strcmp(token[2], "Seq:") == 0) && (strcmp(token[3], "RES:") == 0)) {
                //structureString.assign(token[1], strlen(token[1]));
                token[0] = strtok_r(token[1], "  ", &saveptr1);
                token[2] = strtok_r(saveptr1, "  ", &saveptr2);
                strcpy(output_structure, token[2]);       

                token[2] = strtok_r(saveptr2, "\n", &saveptr2);
                energyString.assign(token[2], strlen(token[2]));
                *output_energy = stod(energyString);
            } else {
                write_log_file(token[0], file, 'E');
                write_log_file(token[1], file, 'E');
                return false;
            }

			/*// Find the structure
			if(boost::regex_search(structureString, match, structureRegex)) {
				structureString = match[1];
			} else {
				write_log_file(token[0], programPath, 'E');
				return false;
			}
			
			token[0] = strtok_r(token[1], "  ", &saveptr1);
			token[2] = strtok_r(saveptr1, "  ", &saveptr2); 
			token[2] = strtok_r(saveptr2, "\n", &saveptr2); 
			
			energyString.assign(token[2], strlen(token[2]));*/

			//Print result
			//std::cout << "Result: " << structureString << " " << energyString << '\n' << std::flush;
	
			// Copy strings
			//strcpy(output_structure, structureString.c_str());
			//*output_energy = stod(energyString);
			//write_log_file("hfold struct", output_structure, 'E');
			return true;
	}
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
	std::string regexResult(boost::regex_replace(structureString, replacedRegex, replacementText));
	strcpy(replaced_structure, regexResult.c_str());
}

void replace_simfold_structure_with_original (char *replaced_structure, char *simfold_structure, int begin, int end) {
	replacementText = "_";
	replacedRegex = "\\.";

	for (int i = begin; i <= end; i++) {
		replaced_structure[i] = simfold_structure[i];
	}

	structureString.assign(replaced_structure, strlen(replaced_structure));
	std::string regexResult(boost::regex_replace(structureString, replacedRegex, replacementText));
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
	sequenceString = boost::regex_replace(sequenceString, sequenceRegex, replacementText);
	std::cout << "Sequence:  " << sequenceString << " Length: " << sequenceString.length() << '\n' << std::flush;

	// Format structure
	structureString = boost::regex_replace(structureString, structureRegex, replacementText);
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

bool get_sequence_structure_interacting (char *fileName, char *sequence_one, char *structure_one, char *sequence_two, char *structure_two) {
	FILE *ioFile;
	size_t len = 0;
	ssize_t read = 0;

	char *sequenceOneBuffer = NULL;
	char *structureOneBuffer = NULL;
	char *sequenceTwoBuffer = NULL;
	char *structureTwoBuffer = NULL;

	std::string sequenceOneString;
	std::string structureOneString;
	std::string sequenceTwoString;
	std::string structureTwoString;

	replacementText = "";
	sequenceRegex = "(\\s|\n|\r)";
	structureRegex = "(\\s|\n|\r|\Z)";

	// Get structure and sequence from file
	ioFile = fopen(fileName, "r");
	if (ioFile != NULL) {
		if((read = getline(&sequenceOneBuffer, &len, ioFile)) == -1) {
			write_log_file("Could not read the first sequence", fileName, 'E');
			return false;
		}
		sequenceOneString.assign(sequenceOneBuffer, read);		

		if((read = getline(&structureOneBuffer, &len, ioFile)) == -1) {
			write_log_file("Could not read the first structure", fileName, 'E');
			return false;
		}
		structureOneString.assign(structureOneBuffer, read);

		if((read = getline(&sequenceTwoBuffer, &len, ioFile)) == -1) {
			write_log_file("Could not read the second sequence", fileName, 'E');
			return false;
		}
		sequenceTwoString.assign(sequenceTwoBuffer, read);		

		if((read = getline(&structureTwoBuffer, &len, ioFile)) == -1) {
			write_log_file("Could not read the second structure", fileName, 'E');
			return false;
		}
		structureTwoString.assign(structureTwoBuffer, read);

	} else {
		write_log_file("Could not open file", fileName, 'E');
		return false;
	}

	// Format sequence
	sequenceOneString = boost::regex_replace(sequenceOneString, sequenceRegex, replacementText);
	sequenceTwoString = boost::regex_replace(sequenceTwoString, sequenceRegex, replacementText);
	std::cout << "Sequence One:  " << sequenceOneString << " Length: " << sequenceOneString.length() << '\n' << std::flush;
	std::cout << "Sequence Two:  " << sequenceTwoString << " Length: " << sequenceTwoString.length() << '\n' << std::flush;

	// Format structure
	structureOneString = boost::regex_replace(structureOneString, structureRegex, replacementText);
	structureTwoString = boost::regex_replace(structureTwoString, structureRegex, replacementText);
	std::cout << "Structure One: " << structureOneString << " Length: " << structureOneString.length() << '\n' << std::flush;
	std::cout << "Structure Two: " << structureTwoString << " Length: " << structureTwoString.length() << '\n' << std::flush;

	// Copy strings
	strcpy(sequence_one, sequenceOneString.c_str());
	strcpy(sequence_two, sequenceTwoString.c_str());
	strcpy(structure_one, structureOneString.c_str());
	strcpy(structure_two, structureTwoString.c_str());

	// Clean up
	fclose(ioFile);
	if (sequenceOneBuffer)
        free(sequenceOneBuffer);
	if (structureOneBuffer)
        free(structureOneBuffer);
	if (sequenceTwoBuffer)
        free(sequenceTwoBuffer);
	if (structureTwoBuffer)
        free(structureTwoBuffer);

	if ((sequenceOneString.length() != structureOneString.length()) || (sequenceTwoString.length() != structureTwoString.length())) {
		write_log_file("Sequence and structure length do not match", fileName, 'E');
		return false;
	}

	return true;
}

bool save_file (char *fileName, char *outputPath, char *sequence, char *restricted, char *structure, double energy, int chosen_method) {
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

void write_log_file(char *message, char *fileName, char option) {
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

void find_sub_sequence_structure (char *input_sequence, char *input_structure, char *output_sequence, char *output_structure, int *begin, int *end) {
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

	if(!boost::regex_search(structureString, match, structureRegex)) {
		return false;
	}

	//substitute all ".()" to "_"
	replacementText = "_";
	structureRegex = "(\\.|\\(|\\))";
	structureString = boost::regex_replace(structureString, structureRegex, replacementText);

	// substitute all [ to (
	replacementText = "(";
	structureRegex = "\\[";
	structureString = boost::regex_replace(structureString, structureRegex, replacementText);

	// substitute all ] to )
	replacementText = ")";
	structureRegex = "\\]";
	structureString = boost::regex_replace(structureString, structureRegex, replacementText);

	strcpy(output_structure, structureString.c_str());
	//std::cout << "New Structure String: " << structureString << '\n' << std::flush;
	return true;
}

int find_no_base_pairs (char *structure) {
	structureString.assign(structure, strlen(structure));
	return std::count(structureString.begin(), structureString.end(), '(');
}

int find_files(char **files_found, char *path) {
	DIR *dir;
	struct dirent *entry;
	int numberOfFiles = 0;
	std::string fileName;
	boost::regex fileNameRegex("(\\w+\\.txt)");

	if ((dir = opendir (path)) != NULL) {
		while ((entry = readdir (dir)) != NULL) {
			fileName = entry->d_name;
			
			if(boost::regex_search(fileName, match, fileNameRegex)) {
				strcpy(files_found[numberOfFiles++], fileName.c_str());
			}
		}
		closedir (dir);

	} else {
		write_log_file("Could not open directory", path, 'E');
		giveup("Could not open directory", path);
	}

	return numberOfFiles;
}

std::string exec(const char* cmd) {
    char buffer[128];
    std::string result = "";
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe)
		return "";
    while (!feof(pipe.get())) {
        if (fgets(buffer, 128, pipe.get()) != NULL)
            result += buffer;
    }
    return result;
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

	call_simfold(SIMFOLD, sequence, NULL, simfold_structure, &simfold_energy, true);
	input_structure_basepair_count = find_no_base_pairs(structure);
	simfold_structure_basepair_count = find_no_base_pairs(simfold_structure);
	ratio = (double) input_structure_basepair_count/simfold_structure_basepair_count;

	//std::cout << "Sparsity ratio is: " << ratio << '\n' << std::flush;

	return ratio;
}

// This function combines two files of both structure and sequence
void file_combiner() {
	int file_number = 1;
	int file_one_length;
	int file_two_length;
	FILE *output_file;
	FILE *file_one;
	FILE *file_two;
	int files_length_one, files_length_two;
	char info[1000];
	char *token[2];
	char *output_path = "\0";
	char *path_one = "./it_test/exons/journal_pone_exons/exon_53_comb/";
	char *path_two = "./it_test/exons/journal_pone_oligos/oligo-exon_53_comb/";
	char *file_one_buffer = "\0";
	char *file_two_buffer = "\0";
	char **files_found_one;
	char **files_found_two;
	
	logFile = fopen(LOGFILEPATH, "w");
	if (logFile == NULL) {
		giveup("Could not open logfile", LOGFILEPATH);
	}

	output_path = (char*) malloc(MAXSLEN*sizeof(char));
	file_one_buffer = (char*) malloc(MAXSLEN*sizeof(char));
	file_two_buffer = (char*) malloc(MAXSLEN*sizeof(char));
	files_found_one = (char**) malloc(MAXSLEN*sizeof(char*));
	files_found_two = (char**) malloc(MAXSLEN*sizeof(char*));
	for (int i=0; i < MAXSLEN; i++) {
		files_found_one[i] = (char*) malloc(sizeof(char) * 255);
		files_found_two[i] = (char*) malloc(sizeof(char) * 255);
	}

	files_length_one = find_files(files_found_one, path_one);
	files_length_two = find_files(files_found_two, path_two);

	write_log_file("Starting Combining", "", 'I');

	for (int i = 0; i < files_length_one; i++) {
		for (int j = 0; j < files_length_two; j++) {
			strcpy(file_one_buffer, path_one);
			strcat(file_one_buffer, files_found_one[i]);
			strcpy(file_two_buffer, path_two);
			strcat(file_two_buffer, files_found_two[j]);

			file_one = fopen(file_one_buffer, "r");
			file_two = fopen(file_two_buffer, "r");

			if (file_one != NULL && file_two != NULL) {
				strcpy(file_two_buffer, files_found_two[j]);
				token[0] = strtok(file_two_buffer, "-");
				token[1] = strtok(NULL, ".");
				
				sprintf(output_path, "./it_test/exons/journal_pone_combined/exon_53/%s-%i-%s.txt", token[0], i+1, token[1], file_number++);
				output_file = fopen(output_path, "w");		

				fseek(file_one, 0, SEEK_END);
   				file_one_length = ftell(file_one);
				fseek(file_one, 0, SEEK_SET);
				
				fseek(file_two, 0, SEEK_END);
   				file_two_length = ftell(file_two);
				fseek(file_two, 0, SEEK_SET);

				fread(file_one_buffer, 1, file_one_length, file_one);
				fread(file_two_buffer, 1, file_two_length, file_two);

				file_one_buffer[file_one_length] = '\0';
				file_two_buffer[file_two_length] = '\0';
				
				fprintf(output_file, file_one_buffer);
				fprintf(output_file, file_two_buffer);

				if (output_file)
					fclose(output_file);
			} else {
				write_log_file("One of the files doesnt exist", "", 'I');
			}

			if (file_one)
				fclose(file_one);
			if (file_two)
				fclose(file_two);
		}

		file_number = 1;
	}

	write_log_file("Ending Combining", "", 'I');

	// Clean up
	for (int i=0; i < MAXSLEN; i++) {
		if (files_found_one[i])
			free(files_found_one[i]);
		if (files_found_two[i])
			free(files_found_two[i]);
	}
	free(output_path);
	free(files_found_one);
	free(files_found_two);
	free(file_one_buffer);
	free(file_two_buffer);
	fclose(logFile);
}

// This function checks text files in two directories and compares them to see if they
// have the exact same contents. Results are recorded in a logfile.
void validator() {
	int line_number;
	int files_length_one, files_length_two;
	double energy_one;
	double energy_two;
	double energy_difference;
	size_t len_one = 0;
	size_t len_two = 0;
	ssize_t read_one = 0;
	ssize_t read_two = 0;
	FILE *file_one;
	FILE *file_two;
	bool file_found;
	char info[1000];
	char *path_one = "./it_test/exons/exon_6_mutated_flank_output_PMO/";
	//char *path_two = "./it_test/exons/exon_6_mutated_flank_output_PMO_temp/";
	char *path_two = "./it_test/exons/exon_6_test/";
	char *file_one_buffer = "\0";
	char *file_two_buffer = "\0";
	char **files_found_one;
	char **files_found_two;
	
	logFile = fopen(LOGFILEPATH, "w");
	if (logFile == NULL) {
		giveup("Could not open logfile", LOGFILEPATH);
	}

	file_one_buffer = (char*) malloc(MAXSLEN*sizeof(char));
	file_two_buffer = (char*) malloc(MAXSLEN*sizeof(char));
	files_found_one = (char**) malloc(MAXSLEN*sizeof(char*));
	files_found_two = (char**) malloc(MAXSLEN*sizeof(char*));
	for (int i=0; i < MAXSLEN; i++) {
		files_found_one[i] = (char*) malloc(sizeof(char) * INF);
		files_found_two[i] = (char*) malloc(sizeof(char) * INF);
	}

	files_length_one = find_files(files_found_one, path_one);
	files_length_two = find_files(files_found_two, path_two);

	write_log_file("Starting Test", "", 'I');

	for (int i = 0; i < files_length_one; i++) {
		file_found = false;

		for (int j = 0; j < files_length_two; j++) {
			if (strcmp(files_found_one[i], files_found_two[j]) == 0) {
				file_found = true;

				strcpy(file_one_buffer, path_one);
				strcat(file_one_buffer, files_found_one[i]);
				strcat(file_one_buffer, "\0");
				strcpy(file_two_buffer, path_two);
				strcat(file_two_buffer, files_found_two[j]);
				strcat(file_two_buffer, "\0");

				file_one = fopen(file_one_buffer, "r");
				file_two = fopen(file_two_buffer, "r");

				line_number = 1; 

				if (file_one != NULL && file_two != NULL) {
					while (((read_one = getline(&file_one_buffer, &len_one, file_one)) != -1) && ((read_two = getline(&file_two_buffer, &len_two, file_two)) != -1)) {

						if (std::strcmp(file_one_buffer, file_two_buffer) == 0) {	
							//sprintf(info, "Line %i matches in file %s", line_number++, files_found_one[i]);
							//write_log_file(info, "", 'I');
						} else {
							//sprintf(info, "Line %d does not match [%s|%s]", line_number++, file_one_buffer, file_two_buffer);
							if (line_number == 4) {
								sscanf(file_one_buffer, "energy: %lf", &energy_one);
								sscanf(file_two_buffer, "energy: %lf", &energy_two);

								energy_difference = std::abs(energy_one) - std::abs(energy_two);
								if (energy_difference < 1 && energy_difference > -1) {
									continue;
								}
							}

							write_log_file("", files_found_one[i], 'E');
							break;
						}
				
						line_number++;
					}
				}

				if (file_one)
					fclose(file_one);
				if (file_two)
					fclose(file_two);
				break;
			}
		}
		
		if (!file_found) {
			write_log_file("file could not be found", files_found_one[i], 'E');
		}
	}

	write_log_file("Ending Test", "", 'I');

	// Clean up
	for (int i=0; i < MAXSLEN; i++) {
		if (files_found_one[i])
			free(files_found_one[i]);
		if (files_found_two[i])
			free(files_found_two[i]);
	}
	free(files_found_one);
	free(files_found_two);
	free(file_one_buffer);
	free(file_two_buffer);
	fclose(logFile);
}

struct Output_Data {
	std::string sequence_one;
	std::string sequence_two;
	std::string final_structure;
	std::string file;
	double energy;

	bool operator > (const Output_Data& structure) const {
        return (energy > structure.energy);
    }
};

struct less_than_energy
{
    inline bool operator() (const Output_Data& struct1, const Output_Data& struct2)
    {
        return (std::abs(struct1.energy) > std::abs(struct2.energy));
    }
};

void compare_output() {
	int line_number;
	size_t len = 0;
	ssize_t read = 0;
	FILE *file_input;
	FILE *file_output;
	FILE *file_map;
	int files_length;
	std::string oligo;
	std::string energy;
	std::string structure;
	std::string sequence_one;
	std::string sequence_two;
	std::string sequence_combined;
	char *path_input = "./it_test/exons/exon_6_mutated_flank_output_PMO/";
	char *path_output = "./it_test/exons/exon_6_mutated_flank_output_PMO/results.txt";
	char *path_map = "./it_test/exons/junk/exon_map_PMO/";
	char *file_buffer = "\0";
	char **files_found;

	// Setup Regex
	static boost::smatch fileMatch;
	static boost::regex fileEnergyRegex("energy: ([^\\n+]*)");
	static boost::regex fileSequenceRegex("(\\w+)XXXXX(\\w+)\n");
	//static boost::regex fileSequenceRegex("(\\w+)\n");
	static boost::regex oligoRegex("exon_\\d+-(\\w+_\\d+).txt");

	Output_Data out_data;

	std::vector<Output_Data> data;
	std::map<std::string, std::string>::iterator it;
	std::map<std::string, std::string> exonMap = {{"GAUUUAUUGGAUCAUUCGUGUACAUCAGGAAGUGGCUCUGGUCUUCCUUUUCUGGUACAAAGAACAGUGGCUCGCCAGAUUACACUGUUGGAGUGUGUCG", "Exon 6"}, {"AACCCAGAUUGCGCCAUUGAACUGCAGCCUGGGCAACAAGAGUGAAACUCCGUCUCAAAAAAAAUAAAAAAUACAAAACAAUAAAACCAGUCCUUCUUCCUUCUUCCAGAGGAGCUUUACAUGUACACUAACAGGCCACGUGUCCCGGAUUGCUGCCCUUCAUGUGAGUUACAAUGUCAUGCAUGCUAAUACUCCAAAGUGGGAGCUAUAUUGCUCAAUCGUUUCUUUUCCCCCUUGUCUUAAACCACAGGAUUUAUUGGAUCAUUCGUGUACAUCAGGAAGUGGCUCUGGUCUUCCUUUUCUGGUACAAAGAACAGUGGCUCGCCAGAUUACACUGUUGGAGUGUGUCGGUAAUUCUUUUUUUUCCUUUCUUUGUGGGUAAUAUGCAAUGUUAGUUUUGUUUUUGAAGUAAAAGAUGGAACUUGGAAAAUCUGCUUUUUGUUACCUGUUGUUUCCUUAUUAAACUUAUGAGAUGUGGUCUCUAUUCUGUAUGACUAAUGUUGGAGCACCUGAACUUGUUUAAGCUCUUAUCAAAGGAGUCCAUCACUAAUGAGUACAGAGAAUGCAGAGUCAAAAGCUUUCAGUUGUAAGUGGACCCUU", "Exon 6 Flank"}, {"GAUUUAUUGGAUCAUUCGUGUACAUCAGGAAGUGGCUCUGGUCUUCCUUUUCUGGUACAAAGAACAGUGGCUCACCAGAUUACACUGUUGGAGUGUGUCG", "Exon 6 Mutated"}, {"AACCCAGAUUGCGCCAUUGAACUGCAGCCUGGGCAACAAGAGUGAAACUCCGUCUCAAAAAAAAUAAAAAAUACAAAACAAUAAAACCAGUCCUUCUUCCUUCUUCCAGAGGAGCUUUACAUGUACACUAACAGGCCACGUGUCCCGGAUUGCUGCCCUUCAUGUGAGUUACAAUGUCAUGCAUGCUAAUACUCCAAAGUGGGAGCUAUAUUGCUCAAUCGUUUCUUUUCCCCCUUGUCUUAAACCACAGGAUUUAUUGGAUCAUUCGUGUACAUCAGGAAGUGGCUCUGGUCUUCCUUUUCUGGUACAAAGAACAGUGGCUCACCAGAUUACACUGUUGGAGUGUGUCGGUAAUUCUUUUUUUUCCUUUCUUUGUGGGUAAUAUGCAAUGUUAGUUUUGUUUUUGAAGUAAAAGAUGGAACUUGGAAAAUCUGCUUUUUGUUACCUGUUGUUUCCUUAUUAAACUUAUGAGAUGUGGUCUCUAUUCUGUAUGACUAAUGUUGGAGCACCUGAACUUGUUUAAGCUCUUAUCAAAGGAGUCCAUCACUAAUGAGUACAGAGAAUGCAGAGUCAAAAGCUUUCAGUUGUAAGUGGACCCUU", "Exon 6 Mutated Flank"}}; //<sequence, oligo/exon>

	// Setup files
	logFile = fopen(LOGFILEPATH, "w");
	if (logFile == NULL) {
		giveup("Could not open logfile", LOGFILEPATH);
	}

	file_buffer = (char*) malloc(MAXSLEN*sizeof(char));
	files_found = (char**) malloc(MAXSLEN*sizeof(char*));
	for (int i=0; i < MAXSLEN; i++) {
		files_found[i] = (char*) malloc(sizeof(char) * 255);
	}

	//Get exon map
	files_length = find_files(files_found, path_map);
	for (int i = 0; i < files_length; i++) {
		oligo.assign(files_found[i], strlen(files_found[i]));
		strcpy(file_buffer, path_map);
		strcat(file_buffer, files_found[i]);
		file_map = fopen(file_buffer, "r");
		
		if ((file_map != NULL) && ((read = getline(&file_buffer, &len, file_map)) != -1)) {
			if(boost::regex_search(oligo, fileMatch, oligoRegex)) {
				oligo = fileMatch[1];
			} else {
				write_log_file("problem with oligo regex", files_found[i], 'E');
				continue;
			}

			file_buffer[strlen(file_buffer)-1] = '\0';		
			exonMap.insert(std::pair<std::string, std::string>(file_buffer, oligo));

			if (file_map)
				fclose(file_map);
		} else {
			write_log_file("Could not open file", file_buffer, 'E');
		}
	}

	files_length = find_files(files_found, path_input);
	file_output = fopen(path_output, "w");	
	write_log_file("Starting Comparison", "", 'I');

	for (int i = 0; i < files_length; i++) {
		energy = "";
		structure = "";
		sequence_one = "";
		sequence_two = "";
		sequence_combined = "";

		strcpy(file_buffer, path_input);
		strcat(file_buffer, files_found[i]);
		file_input = fopen(file_buffer, "r");

		if (file_input != NULL) {
			line_number = 1;

			// Read structure and energy from the output file
			while ((read = getline(&file_buffer, &len, file_input)) != -1) {
				if (line_number == 1) {	
					sequence_combined.assign(file_buffer, strlen(file_buffer));
				} else if (line_number == 3) {	
					structure.assign(file_buffer, strlen(file_buffer));
				} else if (line_number == 4) {	
					energy.assign(file_buffer, strlen(file_buffer));
				} 

				line_number++;
			}

			if (structure.compare("\n") == 0) {
				write_log_file("no structure", files_found[i], 'E');
				continue;
			}			

			// extract information from the structure and energy
			if (line_number == 6) {
				// Find sequences
				if(boost::regex_search(sequence_combined, fileMatch, fileSequenceRegex)) {
					sequence_one = fileMatch[1];
					sequence_two = fileMatch[2];
					//sequence_two = "";
				} else {
					write_log_file("problem with structure regex", files_found[i], 'E');
					continue;
				}
				
				// Find energy
				if(boost::regex_search(energy, fileMatch, fileEnergyRegex)) {
					energy = fileMatch[1];
				} else {
					write_log_file("problem with energy regex", files_found[i], 'E');
					continue;
				}

				out_data = Output_Data();

				// Determine sequence one and two
				it = exonMap.find(sequence_one);
  				if (it != exonMap.end()) {
					out_data.sequence_one = it->second;
				} else {
					write_log_file("Could not find sequence one", files_found[i], 'E');
					continue;
				}

				it = exonMap.find(sequence_two);
  				if (it != exonMap.end()) {
					out_data.sequence_two = it->second;
				} else {
					write_log_file("Could not find sequence two", files_found[i], 'E');
					continue;
				}

				out_data.energy = stod(energy);
				out_data.file = files_found[i];
				structure.replace(sequence_one.length(), linker_length, "XXXXX");
				out_data.final_structure = structure;

				data.push_back(out_data);

			} else {
				write_log_file("Incorrect number of lines", files_found[i], 'E');
			}

			if (file_input)
				fclose(file_input);
			
		} else {
			write_log_file("Could not open file", file_buffer, 'E');
		}
	}

	std::sort(data.begin(), data.end(), less_than_energy());

	// Print all of the data
	/*if (file_output != NULL) {
		for (std::vector<Output_Data>::iterator iter = data.begin(); iter != data.end(); ++iter) {
    		//std::cout << iter->energy << '\n';
			fprintf(file_output, "Energy: %.2f  Sequence 1: %s  Sequence 2: %s  File: %s  Final Structure: %s\n", iter->energy, iter->sequence_one.c_str(), iter->sequence_two.c_str(), iter->file.c_str(), iter->final_structure.c_str());
		}
	} else {
		write_log_file("Cannot open file", path_output, 'E');
	}*/

	// Print the three lowest data
	if (file_output != NULL) {
		for (int i = 0; i < 3; i++) {
    		//std::cout << iter->energy << '\n';
			fprintf(file_output, "Energy: %.2f  Sequence 1: %s  Sequence 2: %s  File: %s  Final Structure: %s\n", data.at(i).energy, data.at(i).sequence_one.c_str(), data.at(i).sequence_two.c_str(), data.at(i).file.c_str(), data.at(i).final_structure.c_str());
		}
	} else {
		write_log_file("Cannot open file", path_output, 'E');
	}

	write_log_file("Ending Comparison", "", 'I');

	// Clean up
	for (int i=0; i < MAXSLEN; i++) {
		if (files_found[i])
			free(files_found[i]);
	}
	free(files_found);
	free(file_buffer);
	fclose(logFile);
	fclose(file_output);
}

void create_oligo_files() {
	int file_number = 1;
	size_t len = 0;
	ssize_t read = 0;
	FILE *output_file;
	FILE *input_file;
	char *input_buffer;
	char *output_buffer = "\0";
	char *input_file_path = "./it_test/exons/journal_pone_oligos/exon_51-oligos.txt";
	char *path_one = "./it_test/exons/journal_pone_oligos/oligo-exon_51/";

	logFile = fopen(LOGFILEPATH, "w");
	if (logFile == NULL) {
		giveup("Could not open logfile", LOGFILEPATH);
	}

	input_file = fopen(input_file_path, "r");

	input_buffer = (char*) malloc(MAXSLEN*sizeof(char));
	output_buffer = (char*) malloc(MAXSLEN*sizeof(char));

	write_log_file("Starting Creation", "", 'I');

	while ((read = getline(&input_buffer, &len, input_file)) != -1) {
		input_buffer[len-1] = '\0';

		strcpy(output_buffer, path_one);
		strcat(output_buffer, "exon_51-oligo_");	
		strcat(output_buffer, std::to_string(file_number++).c_str());
		strcat(output_buffer, ".txt");

		output_file = fopen(output_buffer, "w");
			 
		fprintf(output_file, input_buffer);
		if (output_file)
			fclose(output_file);
	}

	write_log_file("Ending Creation", "", 'I');

	// Clean up
	free(input_buffer);
	free(output_buffer);
	fclose(logFile);
	fclose(input_file);
}
