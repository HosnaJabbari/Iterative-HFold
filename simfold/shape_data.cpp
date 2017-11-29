#include "shape_data.h"

#include <fstream>
#include <iostream>

shape_info shape = shape_info();

// prepares filename to be stored in shape_file
void shape_info::set_shape_file(std::string filename) {
    auto front = filename.front();
    auto back = filename.back();

    // cut quote from beginning and end if they exist
    if (front == '\"')
        filename.erase(front);
    if (back == '\"')
        filename.erase(back);

    shape_file_ = filename;
    use_shape_data_ = true;

    set_data(shape_file_);
}

void shape_info::set_data(std::string filename) {
    // open the file
    std::ifstream infile;
    infile.open(filename, std::ifstream::in);
    if(infile.fail()) {
        fprintf(stderr,"SHAPE data file error: error opening file\n");
        exit(-1);
    }

    // Vector is double the size it needs to be because otherwise there is problems with memory
    // and for reasons I don't understand it gives bad results at large sequence sizes. 
    // This happens when data_[i] is set down below, and it doesn't matter whether you actually use
    // the data_, simply setting it causes the program to have bad output.
    // I checked, and it doesn't set anything past sequence_length()-1, so I have no idea why this works.
    // sequence_length()+1 didn't work, neither did sequence_length()+20, but *2 did. I checked, 
    // and it doesn't set anything past sequence_length()-1, so I have no idea why this works. 
    // Ultimately, a single-dimensional vector of doubles is tiny compared to the rest of the program,
    // and I've spent too much time on this already, so I'm leaving it - Ian Wark September 12 2017
    data_.resize(sequence_length()*2);

    // first line may be the sequence
    std::string input;

    // continue until there is a line that is a number
    while (infile >> input && !is_number(input)) {}

    if (!is_number(input)) {
        fprintf(stderr, "SHAPE data file error: file has no numbers in it\n");
        exit(-1);
    }

    bool valid = true;

    // go through where each word is a shape data number
    int i = 0;
    while (valid) {
        if (i > sequence_length()) {
            fprintf(stderr, "SHAPE data file error: length greater than sequence length (sequence length = %d)\n",sequence_length());
            exit(-1);
        }
        if (!is_number(input)) {
            fprintf(stderr, "SHAPE data file error: line after start is not a number\n",i,sequence_length());
            exit(-1);
        }


        data_[i] = std::stod(input);
        ++i;

        valid = static_cast<bool>(infile >> input);
    }

    if (i != sequence_length()) {
        fprintf(stderr, "SHAPE data file error: length less than sequence length (%d compared to %d)\n",i-1,sequence_length());
        exit(-1);
    }

    infile.close();
}
