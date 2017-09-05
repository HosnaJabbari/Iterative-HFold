#ifndef HFOLD_H_
#define HFOLD_H_

double hfold(char *sequence, char *restricted, char *structure);
double hfold_pkonly(char *sequence, char *restricted, char *structure);
double hfold_iterative(char *sequence, char *restricted, char *structure); // April 3, 2012

void printUsage();

void method1_calculation (char *sequence, char *structure, char *method1_structure, double *method1_energy);
void method2_calculation (char *sequence, char *structure, char *method2_structure, double *method2_energy);
void method3_calculation (char *sequence, char *structure, char *method3_structure, double *method3_energy);
void method4_calculation (char *sequence, char *structure, char *method4_structure, double *method4_energy, int **result);

bool call_HFold (char *programPath, char *input_sequence, char *input_structure, char *output_structure, double *output_energy);
bool call_simfold (char *programPath, char *input_sequence, char *input_structure, char *output_structure, double *output_energy);

bool get_sequence_structure (char *fileName, char *sequence, char *structure);
bool save_file (const char *fileName, char *outputPath, const char *sequence, char *restricted, char *structure, double energy, int chosen_method);
void write_log_file(const char *message, const char *file, const char option);

void segfault_sigaction(int signal, siginfo_t *si, void *arg);


//30 Aug 2017 kevin and Mahyar
//---------------------------------------the enclosed functions are suppose to be the same as the one in Hfold_interacting, if any changes are made, please change that one too--------------------
void remove_structure_intersection(char* G1, char* G0, char* G_p);
int is_invalid_restriction(char* restricted_structure, char* current_structure);
void obtainRelaxedStems(char* G1, char* G2, char* Gresult);
int paired_structure(int i, int j, int *pair_index, int length);
int is_empty_structure(char* input_structure, char* output_structure);
void find_disjoint_substructure(char* structure, std::vector< std::pair<int,int> > &pair_vector);
//---------------------------------------end of enclosed functions --------------------
//30 Aug 2017 kevin and Mahyar
double method1(char *sequence, char *restricted, char *structure);
double method2(char *sequence, char *restricted, char *structure);
double method3(char *sequence, char *restricted, char *structure);
double method4(char *sequence, char *restricted, char *structure);

#endif /*HFOLD_H_*/

