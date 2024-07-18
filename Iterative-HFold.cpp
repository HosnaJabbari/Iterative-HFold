// Iterative HFold files
#include "hotspot.hh"
#include "Result.hh"
#include "cmdline.hh"
#include "hfold_validation.h"
#include "Iterative-HFold.hh"
//simfold files
#include "s_specific_functions.h"
#include "simfold.h"
#include "externs.h"
#include "h_globals.h"
#include "constants.h"
#include "params.h"
#include "common.h"
// a simple driver for the HFold
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stack>
#include <sys/stat.h>
#include <string>
#include <getopt.h>

// this causes less compilation warnings than the #defines
char HFOLD[8] =						"./HFold";
char HFOLD_PKONLY[15] =				"./HFold_pkonly";
char HFOLD_INTERACTING[20] =	    "./HFold_interacting";
char HFOLD_INTERACTING_PKONLY[27] =	"./HFold_interacting_pkonly";
char SIMFOLD[10] =                  "./simfold";

bool exists (const std::string path) {
  struct stat buffer;   
  return (stat (path.c_str(), &buffer) == 0); 
}

void get_input(std::string file, std::string &sequence, std::string &structure){
	if(!exists(file)){
		std::cout << "Input file does not exist" << std::endl;
		exit(EXIT_FAILURE);
	}
	std::ifstream in(file.c_str());
	std::string str;
	int i = 0;
	while(getline(in,str)){
		if(str[0] == '>') continue;
		if(i==0) sequence = str;
		if(i==1) structure = str;
		++i;
	}
	in.close();
}

//check length and if any characters other than ._()
void validateStructure(std::string sequence, std::string structure){
	if(structure.length() != sequence.length()){
		std::cout << " The length of the sequence and corresponding structure must have the same length" << std::endl;
		exit(EXIT_FAILURE);
	}

	//check if any characters are not ._()
	for(char c : structure) {
		if (!(c == '.' || c == '_' || c == '(' || c == ')')){
			std::cout << "Structure must only contain ._(): " << c << std::endl;
			exit(EXIT_FAILURE);
		}
	}
}

//check if sequence is valid with regular expression
//check length and if any characters other than GCAUT
void validateSequence(std::string sequence){

	if(sequence.length() == 0){
		std::cout << "sequence1 or sequence2 is missing" << std::endl;
		exit(EXIT_FAILURE);
	}
  // return false if any characters other than GCAUT -- future implement check based on type
  for(char c : sequence) {
    if (!(c == 'G' || c == 'C' || c == 'A' || c == 'U' || c == 'T')) {
		std::cout << "Sequence contains character " << c << " that is not G,C,A,U, or T." << std::endl;
		exit(EXIT_FAILURE);
    }
  }
}
/* I am changing iterative HFold so that we run 4 methods and choose the structure with the lowest energy as the winner
1) run HFold_PKonly on input structure, if pked, keep the pk bases as input structure and run HFold on them and get the result
2) run HFold on the original input and get the result
3) get the beginning and end of the given structure, give the subsequence and the substructure to simfold restricted as input, get the simfold structure, then give that as input to HFold_PKonly, get the pked bases as new input to give HFold and get the result
4) run simfold restricted with the given input sequence and structure, and get the simfold structure. Only choose part of the structure that contains the original given input structure, then give that to HFold_PKonly as input, get pked bases and run HFold on them and get the result. */
int main (int argc, char **argv) {

	args_info args_info;

	// get options (call getopt command line parser)
	if (cmdline_parser (argc, argv, &args_info) != 0) {
	exit(1);
	}

	std::string seq;
	if (args_info.inputs_num>0) {
	seq=args_info.inputs[0];
	} else {
		if(!args_info.input_file_given) std::getline(std::cin,seq);
	}

	std::string restricted;
    args_info.input_structure_given ? restricted = input_struct : restricted = "";

	std::string fileI;
    args_info.input_file_given ? fileI = input_file : fileI = "";

	std::string fileO;
    args_info.output_file_given ? fileO = output_file : fileO = "";	

	int number_of_suboptimal_structure = args_info.subopt_given ? subopt : 1;
	bool pk_free = args_info.pk_free_given;

	

	if(fileI != ""){
		
		if(exists(fileI)){
			get_input(fileI,seq,restricted);
		}
		if(seq == ""){
			std::cout << "sequence is missing from file" << std::endl; 
		}
		
	}

	int n = seq.length();

	validateSequence(seq);

	if(restricted != "") validateStructure(seq,restricted);

	// //set up for RNA so we can use this for building hotspot
	// if(!args_info.input_structure_given){
	// 	char config_file[200];
	// 	strcpy (config_file, SIMFOLD_HOME "/params/multirnafold.conf");
	// 	int dna_or_rna = RNA;
	// 	double temperature = 37.0;
	// 	init_data ("./simfold", config_file, dna_or_rna, temperature);
	// 	fill_data_structures_with_new_parameters ( SIMFOLD_HOME "/params/turner_parameters_fm363_constrdangles.txt");
	// 	fill_data_structures_with_new_parameters ( SIMFOLD_HOME "/params/parameters_DP09.txt");
	// }

	char config_file[400];
	strcpy (config_file, SIMFOLD_HOME "/params/multirnafold.conf");

	//what to fold: RNA or DNA
	int dna_or_rna= RNA;
	// represents degrees Celsius
	double temperature = 37.0;
	// initialize the thermodynamic parameters
	init_data ("./simfold", config_file, dna_or_rna, temperature);

	fill_data_structures_with_new_parameters ( SIMFOLD_HOME "/params/turner_parameters_fm363_constrdangles.txt");

	// in HotKnots and ComputeEnergy package the most up-to-date parameters set is DP09.txt
	// so we add it here
	fill_data_structures_with_new_parameters ( SIMFOLD_HOME "/params/parameters_DP09.txt");

	std::vector<Hotspot> hotspot_list;

	char sequence[n+1];
	strcpy(sequence,seq.c_str());
	// Hotspots

	
	if(restricted != ""){
		Hotspot hotspot(0,restricted.length()-1,restricted.length());
		hotspot.set_structure(restricted);
		hotspot_list.push_back(hotspot);
	}
	else {
		get_hotspots(sequence, hotspot_list,number_of_suboptimal_structure);
	}

	// Data structure for holding the output
	std::vector<Result> result_list;

	// Iterate through all hotspots or the single given input structure
	for(int i = 0;i<hotspot_list.size();++i){

		char structure[n+1];
		std::string struc = hotspot_list[i].get_structure();
		strcpy(structure,struc.c_str());

		char final_structure[MAXSLEN];
		double final_energy = INF;
		int method_chosen = -1;
		if(!pk_free){
			char *method1_structure = (char*) malloc(sizeof(char) * MAXSLEN);
			char *method2_structure = (char*) malloc(sizeof(char) * MAXSLEN);
			char *method3_structure = (char*) malloc(sizeof(char) * MAXSLEN);
			char *method4_structure = (char*) malloc(sizeof(char) * MAXSLEN);

			double *method1_energy = (double*) malloc(sizeof(double) * INF);
			double *method2_energy = (double*) malloc(sizeof(double) * INF);
			double *method3_energy = (double*) malloc(sizeof(double) * INF);
			double *method4_energy = (double*) malloc(sizeof(double) * INF);

			*method1_energy = INF;
			*method2_energy = INF;
			*method3_energy = INF;
			*method4_energy = INF;

			method1_structure[0] = '\0';
			method2_structure[0] = '\0';
			method3_structure[0] = '\0';
			method4_structure[0] = '\0';
			final_structure[0] = '\0';

			// printf("method1\n");
			*method1_energy = method1(sequence, structure, method1_structure);
			// printf("method2\n");
			*method2_energy = method2(sequence, structure, method2_structure);
			// printf("method3\n");
			*method3_energy = method3(sequence, structure, method3_structure);
			// printf("method4\n");
			*method4_energy = method4(sequence, structure, method4_structure);

				//We ignore non-negetive energy, only if the energy of the input sequnces are non-positive!
				if (*method1_energy < final_energy) {
				//if (*method1_energy < final_energy && *method1_energy != 0) {
						final_energy = *method1_energy;
						strcpy(final_structure, method1_structure);
						method_chosen = 1;
				}

				if (*method2_energy < final_energy) {
				//if (*method2_energy < final_energy && *method2_energy != 0) {
						final_energy = *method2_energy;
						strcpy(final_structure, method2_structure);
						method_chosen = 2;
				}

				if (*method3_energy < final_energy) {
				//if (*method3_energy < final_energy && *method3_energy != 0) {
						final_energy = *method3_energy;
						strcpy(final_structure, method3_structure);
						method_chosen = 3;
				}

				if (*method4_energy < final_energy) {
				//if (*method4_energy < final_energy && *method4_energy != 0) {
						final_energy = *method4_energy;
						strcpy(final_structure, method4_structure);
						method_chosen = 4;
				}


			if (final_energy == INF || method_chosen == -1) {
				fprintf(stderr, "ERROR: could not find energy\n");
				fprintf(stderr, "SEQ: %s\n",sequence);
				fprintf(stderr, "Structure: %s\n",structure);
				exit(6);
			}

			// Clean up
			free(method1_structure);
			free(method2_structure);
			free(method3_structure);
			free(method4_structure);
			free(method1_energy);
			free(method2_energy);
			free(method3_energy);
			free(method4_energy);
		} 
		else{	
			call_simfold(SIMFOLD, sequence, structure, final_structure, &final_energy);
		}

		std::string final(final_structure);
		Result result(seq,hotspot_list[i].get_structure(),hotspot_list[i].get_energy(),final,final_energy,method_chosen);
		result_list.push_back(result);
	}
	
	Result::Result_comp result_comp;
	std::sort(result_list.begin(), result_list.end(),result_comp );

	int number_of_output = 1;

	if(number_of_suboptimal_structure != 1){
			number_of_output = std::min( (int) result_list.size(),number_of_suboptimal_structure);
	}

	//output to file
	if(fileO != ""){
		std::ofstream out(fileO);
		out << sequence << std::endl;
		for (int i=0; i < number_of_output; i++) {
			out << "Restricted_" << i << ": " << result_list[i].get_restricted() << std::endl;;
			out << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")" << std::endl;
			if(args_info.verbose_given){
				out << "Method: " << result_list[i].get_method_chosen() << std::endl;
			}	
		}

	}else{
		//kevin: june 22 2017
		//Mateo: Sept 13 2023
		//changed format for ouptut to stdout
		std::cout << seq << std::endl;
		if(result_list.size() == 1){
			// std::cout << "Restricted_" << 0 << ": " << result_list[0].get_restricted() << std::endl;;
			std::cout << result_list[0].get_final_structure() << " (" << result_list[0].get_final_energy() << ")" << std::endl;
			
			if(args_info.verbose_given){
				std::cout << "Method: " << result_list[0].get_method_chosen() << std::endl;
			}
		}
		else{
			for (int i=0; i < number_of_output; i++) {
				if(result_list[i].get_final_structure() == result_list[i-1].get_final_structure()) continue;
				std::cout << "Restricted_" << i << ": " << result_list[i].get_restricted() << std::endl;;
				std::cout << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")" << std::endl;
				if(args_info.verbose_given){
					std::cout << "Method: " << result_list[i].get_method_chosen() << std::endl;
				}
			}
		}
	}
	cmdline_parser_free(&args_info);
	return 0;
}



void printUsage(){
	printf("Usage ./Iterative-HFold --s <sequence> --r <structure> [--o </path/to/file>]\n");
	printf("or\n");
	printf("Usage ./Iterative-HFold --i </path/to/file> [--o </path/to/file>]\n");
	printf ("  Restricted structure symbols:\n");
	printf ("    () restricted base pair\n");
	printf ("    _ no restriction\n");
	printf("Example:\n");
	printf("./Iterative-HFold --s \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" --r \"(____________________________)\"\n");
	printf("./Iterative-HFold --i \"/home/username/Desktop/myinputfile.txt\" --o \"/home/username/Desktop/some_folder/outputfile.txt\"\n");
	printf("Please read README for more details\n");
}


void segfault_sigaction(int signal, siginfo_t *si, void *arg) {
	//Need to expand on this.
    fprintf(stderr, "Caught segfault at address %p\n", si->si_addr);
    exit(si->si_errno);
}



//Calling HFold: programPath = HFOLD
//Calling HFold_PKonly: programPath = HFOLD_PKONLY
bool call_HFold (char *programPath, char *input_sequence, char *input_structure, char *output_structure, double *output_energy) {
	char config_file[400];
	strcpy (config_file, SIMFOLD_HOME "/params/multirnafold.conf");

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

	fill_data_structures_with_new_parameters ( SIMFOLD_HOME "/params/turner_parameters_fm363_constrdangles.txt");

	// in HotKnots and ComputeEnergy package the most up-to-date parameters set is DP09.txt
	// so we add it here
	fill_data_structures_with_new_parameters ( SIMFOLD_HOME "/params/parameters_DP09.txt");

        if (strcmp(programPath, HFOLD) == 0) {
            *output_energy = hfold(input_sequence, input_structure, output_structure);
        }
        else if (strcmp(programPath, HFOLD_PKONLY) == 0){
                *output_energy = hfold_pkonly(input_sequence, input_structure, output_structure);
        }else{
                printf("Error: invalid arguments are given: %s \nValid aurgumnets are: HFOLD and HFOLD_PKONLY\n");
                return false;
         }
	if(is_invalid_restriction(input_structure,output_structure)){
		fprintf(stderr,"ERROR!!! There is something wrong with the structure, doesn't match restricted\n");
        fprintf(stderr,"  %s\n  %s\n  %s\t%.2lf\n", input_sequence, input_structure, output_structure, *output_energy);
        fprintf(stderr,"ERROR!!! There is something wrong with the structure, doesn't match restricted\n");
		exit(11);
	}
	return true;

}

bool call_simfold (char *programPath, char *input_sequence, char *input_structure, char *output_structure, double *output_energy) {
        std::string result = "";

		char config_file[400];
		strcpy (config_file, SIMFOLD_HOME "/params/multirnafold.conf");

		double temperature;
		temperature = 37;
		init_data ("./simfold", config_file, RNA, temperature);

        fill_data_structures_with_new_parameters (SIMFOLD_HOME "/params/turner_parameters_fm363_constrdangles.txt");
	// when I fill the structures with DP09 parameters, I get a segmentation fault for 108 base sequence!!!!
	// So I chopped the parameter set to only hold the exact number as the turner_parameters_fm363_constrdangles.txt,
	// but still getting seg fault!
	fill_data_structures_with_new_parameters (SIMFOLD_HOME "/params/parameters_DP09_chopped.txt");

	*output_energy = simfold_restricted (input_sequence, input_structure, output_structure);
//	printf ("Call_Simfold_RES( can be called by different methods): %s  %.2lf\n", output_structure, output_energy);
	if(is_invalid_restriction(input_structure,output_structure)){
			fprintf(stderr,"ERROR!!! There is something wrong with the structure, doesn't match restricted\n");
			fprintf(stderr,"  %s\n  %s\n  %s\t%.2lf\n", input_sequence, input_structure, output_structure, *output_energy);
			fprintf(stderr,"ERROR!!! There is something wrong with the structure, doesn't match restricted\n");
			exit(11);
	}
	return true;
}





//kevin 28 Aug
// Function to remove all spaces from a given string
void removeSpaces(char *str)
{
    // To keep track of non-space character count
    int count = 0;

    // Traverse the given string. If current character
    // is not space, then place it at index 'count++'
    for (int i = 0; str[i]; i++)
        if (str[i] != ' ' || str[i] != '\n' || str[i] != '\t' || str[i] != '\r')
            str[count++] = str[i]; // here count is incremented
    str[count-1] = '\0';
}


//---------------------------------------this function is suppose to be the same as the one in Hfold_interacting, if any changes are made, please change that one too--------------------
//kevin 19 july (verified by Mahyar 19 july 2017)
//count+1 when open bracket, count-1 when close bracket
//whgen count is 0, we have a substructure
//assume valid structure
void find_disjoint_substructure(char* structure, std::vector< std::pair<int,int> > &pair_vector){
	int length = strlen(structure);
	int count = 0;
	int first_time = 1; //flag for first time getting a open bracket for substructure
	int i = 0;
	int j = 0;
	for(int k=0; k<length;k++){
		if(structure[k] == '(' || structure[k] == '['){
			if(first_time && count == 0){
				first_time = 0;
				i = k;
			}
			count += 1;

		}else if(structure[k] == ')' || structure[k] == ']'){
			count -= 1;
			j = k;
			if(count == 0){
				std::pair <int,int> ij_pair (i,j);
				pair_vector.push_back(ij_pair);
				first_time = 1;
			}
		}
	}
}

//---------------------------------------this function is suppose to be the same as the one in Hfold_interacting, if any changes are made, please change that one too--------------------
//31 Aug 2017 kevin and Mahyar
//input[i] is _ or . and it did not turn into a . in the output structure, then it is not empty
int is_empty_structure(char* input_structure, char* output_structure){
	for(int i=0; i<strlen(input_structure);i++){
		if(input_structure[i] != output_structure[i]){
			if((input_structure[i] == '_' || input_structure[i] == '.' ) && output_structure[i] != '.'){
				return 0;
			}
		}
	}
	return 1;
}


//---------------------------------------this function is suppose to be the same as the one in Hfold_interacting, if any changes are made, please change that one too--------------------
//kevin 18 July
int paired_structure(int i, int j, int *pair_index, int length){
	if(i >= 0 && j < length && (pair_index[i] == j) && (pair_index[j] == i) ){
		return 1;
	}
	return 0;
}

//---------------------------------------this function is suppose to be the same as the one in Hfold_interacting, if any changes are made, please change that one too--------------------
//kevin 18 July
void obtainRelaxedStems(char* G1, char* G2, char* Gresult){
	int length = strlen(G1);
	int G1_pair[length];
	int G2_pair[length];

	//Gresult <- G1
	strncpy(Gresult,G1,length);

	detect_original_pairs(G1, G1_pair);
	detect_original_pairs(G2, G2_pair);

	//for(int d=0;d<length;d++){
	//	printf("%c %d\n",G2[d],G2_pair[d]);
	//}

	int i = 0;
	int j = 0;

	for(int k=0;k<length;++k){
		if(G2_pair[k] > -1){
			i = k;
			j = G2_pair[k];
			if(i < j){ //for each ij in G2
				if( (G1[i] != G2[i]) && (G1[j] != G2[j]) ){//if ij not in G1
					//include stacking base pairs
					if(paired_structure(i-1,j+1,G1_pair,length)){
						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include bulges of size 1
					}else if(paired_structure(i-2,j+1,G1_pair,length) || paired_structure(i-1,j+2,G1_pair,length)){
						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include loops of size 1x1
					}else if(paired_structure(i-2,j+2,G1_pair,length)){
						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include loops of size 1x2 or 2x1
					}else if( paired_structure(i-3,j+2,G1_pair,length) || paired_structure(i-2,j+3,G1_pair,length)){
						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					}
				}
			}
		}
	}

	for(int k=length-1;k>=0;--k){
		if(G2_pair[k] > -1){
			i = k;
			j = G2_pair[k];
			if(i < j){ //for each ij in G2
				if( (G1[i] != G2[i]) && (G1[j] != G2[j]) ){//if ij not in G1
					//include stacking base pairs
					if(paired_structure(i+1,j-1,G1_pair,length) ){
						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
						
					//include bulges of size 1
					}else if(paired_structure(i+1,j-2,G1_pair,length) || paired_structure(i+2,j-1,G1_pair,length) ){
						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include loops of size 1x1
					}else if(paired_structure(i+2,j-2,G1_pair,length) ){
						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					//include loops of size 1x2 or 2x1
					}else if(paired_structure(i+2,j-3,G1_pair,length) || paired_structure(i+3,j-2,G1_pair,length) ){
						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
						G1_pair[i] = j;
						G1_pair[j] = i;
					}
				}
			}
		}
	}
}




//---------------------------------------this function is suppose to be the same as the one in Hfold_interacting, if any changes are made, please change that one too--------------------
//kevin 30 Aug 2017
//check if the computed structure matches the restricted structure
int is_invalid_restriction(char* restricted_structure, char* current_structure){
	std::string openBracketArray ("({[");
	std::string closeBracketArray (")}]");

	for (int i=0; i < strlen(restricted_structure); i++){
        if(restricted_structure[i] != '_' && restricted_structure[i] != current_structure[i]){
			if( (openBracketArray.find_first_of(restricted_structure[i]) != -1) && ((openBracketArray.find_first_of(current_structure[i]) != -1)) ){
				continue;
			}else if ( (closeBracketArray.find_first_of(restricted_structure[i]) != -1) && ((closeBracketArray.find_first_of(current_structure[i]) != -1)) ){
				continue;
			}else{
				return 1;
			}
		}

    }
	return 0;
}


//30 Aug 2017 kevin and Mahyar
double method1(char *sequence, char *restricted, char *structure){
	//printf("inside method1 restricrted:\n%s\n",restricted);
	double energy = 0;
	call_HFold(HFOLD, sequence, restricted, structure, &energy);
	//printf("method1: %s\n",structure);
	return energy;
}


//30 Aug 2017 kevin and Mahyar
double method2(char *sequence, char *restricted, char *structure){
	//printf("inside method2 restricrted:\n%s\n",restricted);
	double energy = 0;

	call_HFold(HFOLD_PKONLY, sequence, restricted, structure, &energy);
	//printf("inside method2 structure:\n%s\n",structure);
	//printf("restricted: %s\nstructure: %s\n",restricted,structure);
	if(is_empty_structure(restricted,structure)){
		//printf("is empty\n");
		return energy;
	}else{
		//printf("is not empty\n");
		char G_prime[strlen(structure)];
		remove_structure_intersection(structure,restricted, G_prime);
		//printf("G_prime: %s\n",G_prime);
		energy = method1(sequence, G_prime, structure);
		return energy;
	}

}

//30 Aug 2017 kevin and Mahyar
double method3(char *sequence, char *restricted, char *structure){
	double energy = 0;
	int length = strlen(sequence);
	char simfold_structure[length];

	call_simfold(SIMFOLD, sequence, restricted, simfold_structure, &energy);

	//^ G' simfold_structure <- SimFold(S sequence, G restricted)
	char G_updated[length+1];
	G_updated[length] = '\0';
	obtainRelaxedStems(restricted ,simfold_structure, G_updated);
	//^Gupdated G_updated<- ObtainRelaxedStems(G restricted,G' simfold_structure)
	energy = method2(sequence, G_updated, structure);

	return energy;
}

//30 Aug 2017 kevin and Mahyar
double method4(char *sequence, char *restricted, char *structure){
	int KEVIN_DEBUG = 0;
	double energy = 0;
	int length = strlen(sequence);
	char G_updated[length+1];
	int k = 1;
	//^k <- 1
	strcpy(G_updated, restricted);
	//^Gupdated <- G
	std::vector< std::pair<int,int> > disjoint_substructure_index; //contain pair(i,j) in a vector, ij are index of begin and end of substructure
	find_disjoint_substructure(restricted, disjoint_substructure_index);
	//^get disjoint substructure
	int i = 0;
	int j = 0;

	//31 Aug 2017 kevin and Mahyar
	if(disjoint_substructure_index.size() == 0){
		return INF;
	}

	for(auto current_substructure_index : disjoint_substructure_index){
		i = current_substructure_index.first;
		j = current_substructure_index.second;
		char subsequence[length+1];
		char substructure[length+1];
		char simfold_structure[length+1];

		strncpy(subsequence, sequence+i,j-i+1);
		subsequence[j-i+1] = '\0';
		//^Sk
		strncpy(substructure, restricted+i,j-i+1);
		substructure[j-i+1] = '\0';
		//^Gk

		call_simfold(SIMFOLD, subsequence, substructure, simfold_structure, &energy);

		//^ SimFold(Sk,Gk,Gk')
		char Gp_k_updated[length];
		obtainRelaxedStems(substructure, simfold_structure, Gp_k_updated);
		//^obtainRelaxedStems(Gk,Gk',G'kupdated)

		int m = 0; //index for going through Gp_k_updated
        for(int k =i;k<j;k++){
			if(G_updated[k] != Gp_k_updated[m]){
				G_updated[k] = Gp_k_updated[m];
			}
			m++;
		}
		//^Gupdated <- Gupdated U G'kupdated
	}
	energy = method2(sequence, G_updated, structure);

	return energy;
}

//---------------------------------------this function is suppose to be the same as the one in Hfold_interacting, if any changes are made, please change that one too--------------------
// Aug 31, 2017 kevin and Mahyar
//does G_p = G1-G
void remove_structure_intersection(char* G1, char* G0, char* G_p){
	strcpy(G_p,G1);
	for(int i=0; i< strlen(G1); i++){
		if(G_p[i] == '.'){
			G_p[i] = '_';
		}
		if (G1[i] == G0[i]){
			G_p[i] = '_';
		}
		if (G_p[i] == '[' || G_p[i] == '{'){
			G_p[i] = '(';
		}
		if (G_p[i] == ']' || G_p[i] == '}'){
			G_p[i] = ')';
		}
	}
}
