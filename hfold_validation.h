//created by Kevin Chang 19 June 2017

//change filename.txt to filename_timestamp.txt
void addTimestamp(char** path);

//check if sequence is valid with regular expression
//check length and if any characters other than GCAUT
bool validateSequence(const char *string);

//check if structure is valid
//check if length match with sequence
//check if any characters other than ._(){}[]
//check if pk-free
bool validateStructure(char* structure, char* sequence);

//change all {} [] to () in structure
void replaceBrackets(char* structure);

//if output_path is not a path -> take base path of input file path (if exist) and concat with output_path
//if output_path is a path -> do nothing
void addPath(char** output_path, char* input_path);

//check if input file is in correct format and store value into corresponding variable
//return true on success, false on fail
bool validateInteractingInputFile(char* path, char* seq1, char* struc1, char* seq2, char* struc2);

bool validateHFOLDInputFile(char* path, char* seq1, char* struc1);