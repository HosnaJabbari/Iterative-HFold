#ifndef HFOLD_H_
#define HFOLD_H_

double hfold(char *sequence, char *restricted, char *structure);
double hfold_iterative(char *sequence, char *restricted, char *structure); // April 3, 2012

static void *threadFunction(void *arg);


void printUsage();

void method1_calculation (char *sequence, char *structure, char *method1_structure, double *method1_energy);
void method2_calculation (char *sequence, char *structure, char *method2_structure, double *method2_energy);
void method3_calculation (char *sequence, char *structure, char *method3_structure, double *method3_energy);
void method4_calculation (char *sequence, char *structure, char *method4_structure, double *method4_energy, int **result);

bool call_HFold (char *programPath, char *input_sequence, char *input_structure, char *output_structure, double *output_energy);
bool call_simfold (char *programPath, char *input_sequence, char *input_structure, char *output_structure, double *output_energy, bool reattempt_run);

void replace_simfold_partial_structure_with_original (char *input_structure, char *simfold_structure, char *replaced_structure, int begin, int end);
void replace_simfold_structure_with_original (char *replaced_structure, char *simfold_structure, int begin, int end);
bool get_sequence_structure (char *fileName, char *sequence, char *structure);
bool get_sequence_structure_interacting (char *fileName, char *sequence_one, char *structure_one, char *sequence_two, char *structure_two);
bool save_file (const char *fileName, char *outputPath, const char *sequence, char *restricted, char *structure, double energy, int chosen_method);
void write_log_file(const char *message, const char *file, const char option);

bool find_each_substructure (char *input_structure, int begin, int end, int **output);
void find_independant_structures (char *structure, int begin, int end, int *B, int *Bp);
void find_sub_sequence_structure (const char *input_sequence, char *input_structure, char *output_sequence, char *output_structure, int *begin, int *end);
bool find_new_structure (char *input_structure, char *output_structure);
int find_no_base_pairs (char *structure);

double find_structure_sparsity (char *sequence, char *structure); // Unused. keep for now.
std::string exec(const char* cmd);

void segfault_sigaction(int signal, siginfo_t *si, void *arg);


#endif /*HFOLD_H_*/

