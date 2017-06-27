Description:
These methods are base on hierarchical folding hypothesis

Supported OS: *nix operating system

How to install:
type "make"

How to use:
    
    Arguments:
        HFold:
            -s <sequence>
            -r <structure>
            -i </path/to/file>
            -o </path/to/file>

        Remarks:
            make sure the <arguments> are enclosed in "", for example -r "..().." instead of -r ..()..
            input file for -i must be .txt
            if -i is provided with just a file name without a path, it is assuming the file is in the diretory where the executable is called
            if -o is provided with just a file name without a path, the output file will be generated in the diretory where the executable is called
            if -o is provided with just a file name without a path, and if -i is provided, then the output file will be generated in the directory where the input file is located
    
    Sequence requirements:
        containing only characters GCAUT

    Structure requirements:
        -pseudoknot free
        -containing only characters ._(){}[]
        Remarks:
            Restricted structure symbols:
                () restricted base pair
                _ no restriction


    Input file requirements:
            Line1: Sequence
            Line2: Structure
        sample:
            GCAACGAUGACAUACAUCGCUAGUCGACGC
            (____________________________)

    Example:
    #assume you are in the directory where the HFold/HFold_interacting/HFold_iterative executable is loacted
    ./HFold -i "/home/username/Desktop/myinputfile.txt"
    ./HFold -i "/home/username/Desktop/myinputfile.txt" -o "outputfile.txt"
    ./HFold -i "/home/username/Desktop/myinputfile.txt" -o "/home/username/Desktop/some_folder/outputfile.txt"
    ./HFold -s "GCAACGAUGACAUACAUCGCUAGUCGACGC" -r "(____________________________)"
    ./HFold -s "GCAACGAUGACAUACAUCGCUAGUCGACGC" -r "(____________________________)" -o "outputfile.txt"
    ./HFold -s "GCAACGAUGACAUACAUCGCUAGUCGACGC" -r "(____________________________)" -o "/home/username/Desktop/some_folder/outputfile.txt"

Exit code:
    0       success
    1	    invalid argument error 
    3	    thread error
    4       i/o error
    5       pipe error
    error code with special meaning: http://tldp.org/LDP/abs/html/exitcodes.html
    2	    Misuse of shell builtins (according to Bash documentation)
    126	    Command invoked cannot execute
    127	    "command not found"
    128	    Invalid argument to exit	
    128+n	Fatal error signal "n"
    130	    Script terminated by Control-C
    255	    Exit status out of range (range is 0-255)