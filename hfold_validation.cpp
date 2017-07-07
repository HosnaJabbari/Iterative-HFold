//created by Kevin Chang 19 June 2017
#define _XOPEN_SOURCE_EXTENDED  1
#include <libgen.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <stack>

#include "externs.h"


#include <ctime>
#include <regex>
#include <regex.h>

#include "hfold_validation.h"



static std::regex structureRegex;

//change all {} [] to () in structure
void replaceBrackets(char* structure){
	std::string temp = structure;
	std::string replace = "(";
	structureRegex = "[[{]";
	temp = std::regex_replace(temp, structureRegex, replace);
	replace = ")";
	structureRegex = "[\]}]";
	temp = std::regex_replace(temp, structureRegex, replace);
	//std::cout <<  temp << '\n' << std::flush;
	temp.copy(structure,temp.length(),0);
}

//check if structure is valid
//check if length match with sequence
//check if any characters other than ._(){}[]
//check if pk-free
bool validateStructure(char* structure, char* sequence){
	int    status;
    regex_t    re;

	if(strlen(structure) > MAXSLEN){
		printf("Invalid length\n");
		return false;
	}

	if(strlen(structure) <= 0){
		return false;
	}

	//printf("strlen: %d %d\n",strlen(structure), strlen(sequence));
	if(strlen(structure) != strlen(sequence)){
		printf("Length of sequence and corresponding structure must have same length\n");
		return false;
	}
	
	std::regex expr{"([^\\[\\]._{}()])"};
	//std::cout << std::regex_search(s, expr) << '\n';

    if (std::regex_search(structure, expr) != 0) {
		printf("Structure must only contain ._(){}[] \n");
        return false;      
    }
  
	std::string openBracketArray ("({[");
	std::string closeBracketArray (")}]");
	std::stack<char> mystack;

	for(int i=0; i<strlen(structure);i++){
		//printf("d: %c\n",structure[i]);
		if(openBracketArray.find_first_of(structure[i]) != -1){ //is open bracket
			mystack.push(structure[i]);
		}else if(closeBracketArray.find_first_of(structure[i]) != -1){ //is close bracket
			if(closeBracketArray.find_first_of(structure[i]) != openBracketArray.find_first_of(mystack.top())){ //if current bracket is not corresponding bracket of what we popped
				printf("Brackets must match and is pseudoknot free\n");
				return false;
			}
			mystack.pop();
		}
	}
	if(mystack.empty() == false){
		//printf("fail isempty\n");
		return false;
	}
	return true;      
    
	
}

//check if sequence is valid with regular expression
//check length and if any characters other than GCAUT
bool validateSequence(const char* string){
    int    status;
    regex_t    re;

	if(strlen(string) > MAXSLEN){
		return false;
	}

	if(strlen(string) <= 0){
		return false;
	}

    if (regcomp(&re, "[^GCAUT]", REG_EXTENDED|REG_NOSUB) != 0) {
        return false;     
    }
    status = regexec(&re, string, (size_t) 0, NULL, 0);
    regfree(&re);
    if (status != 0) {
        return true;      
    }
    return false;
}


//if output_path is not a path -> take base path of input file path (if exist) and concat with output_path
//if output_path is a path -> do nothing
void addPath(char** output_path, char* input_path){
	std::string temp_out_path = *output_path;
	std::string temp_in_path = input_path;
	std::size_t out_path_found = temp_out_path.rfind("/");
	
	if(out_path_found == std::string::npos){ //if out path does not contain '/'
		std::size_t in_path_found = temp_in_path.rfind("/");
		if(in_path_found != std::string::npos){ //if in path contain '/'
			std::string base_path = temp_in_path.substr(0,in_path_found+1);
			std::cout << base_path << '\n';
			temp_out_path = base_path + temp_out_path;
			//std::cout << temp_out_path << '\n';
			temp_out_path.copy(*output_path, temp_out_path.length(), 0);
		}
	}
	
}

//assume '.' exist, if not found then do nothings
//change filename.txt to filename_timestamp.txt
void addTimestamp(char** path){
	std::string temp = *path;
	//std::cout << temp << '\n';
	std::size_t found = temp.find_last_of("."); 
  	
	time_t rawtime;
	struct tm * timeinfo;
	char buffer [80];
	time (&rawtime);
	timeinfo = localtime (&rawtime);
	strftime (buffer,80,"_%F_%X",timeinfo);
	//june 22 2017 http://www.cplusplus.com/reference/string/string/find_last_of/
	if(found != std::string::npos){
		temp.insert(found,buffer);
		temp.copy(*path,temp.length(),0);
	}
}

//check if input file is in correct format and store value into corresponding variable
//return true on success, false on fail
bool validateInteractingInputFile(char* path, char* seq1, char* struc1, char* seq2, char* struc2){
    //printf("path: %s\n", path);
    FILE* fp;
    char line[MAXSLEN];

    fp = fopen(path, "r");
    if(!fp){
        printf("File not found\n");
        return false;
    }
    fscanf(fp,"%s\n%s\n%s\n%s",seq1,struc1,seq2,struc2);
    //printf("%s|%s|%s|%s|\n",seq1,struc1,seq2,struc2);
    fclose(fp);

    if(!(strlen(seq1) < strlen(seq2))){
        printf("Line 1 must be shorter than Line 3\n");
        return false;
    }

    if(!(validateSequence(seq1))){
        printf("Line 1 is invalid\n");
        return false;
    }

    if(!(validateSequence(seq2))){
        printf("Line 3 is invalid\n");
        return false;
    }

    if(!(validateStructure(struc1,seq1))){
        printf("Line 2 is invalid\n");
        return false;
    }

    if(!(validateStructure(struc2,seq2))){
        printf("Line 4 is invalid\n");
        return false;
    }

    
    return true;
}

//check if input file is in correct format and store value into corresponding variable
//return true on success, false on fail
bool validateHFOLDInputFile(char* path, char* seq1, char* struc1){
    //printf("path: %s\n", path);
    FILE* fp;
    char line[MAXSLEN];

    fp = fopen(path, "r");
    if(!fp){
        printf("File not found\n");
        return false;
    }
    fscanf(fp,"%s\n%s\n",seq1,struc1);
    
    fclose(fp);

    if(!(validateSequence(seq1))){
        printf("Line 1 is invalid\n");
        return false;
    }

    if(!(validateStructure(struc1,seq1))){
        printf("Line 2 is invalid\n");
        return false;
    }

    return true;
}

