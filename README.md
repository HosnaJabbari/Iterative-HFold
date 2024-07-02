# HFold Iterative

#### Description:
Software implementation of Iterative HFold.      
Iterative HFold is an algorithm for predicting the pseudoknotted secondary structures of RNA using relaxed Hierarchical Folding. 

Paper: Jabbari, H., Condon, A. A fast and robust iterative algorithm for prediction of RNA pseudoknotted secondary structures. BMC Bioinformatics 15, 147 (2014). https://doi.org/10.1186/1471-2105-15-147 (https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-147)

On the dataset tested in this paper, Iterative HFold generally has better accuracy that its predecessor, [HFold](https://github.com/HosnaJabbari/HFold).

#### Supported OS: 
Linux, macOS

### Conda Package:
```
conda install -c uvic-cobra iterative-hfold
```
Works for Linux and macOS

### Source code Installation:  
Requirements: A compiler that supports C++11 standard (tested with g++ version 4.9.0 or higher), Pthreads, and CMake version 3.1 or greater.    

[CMake](https://cmake.org/install/) version 3.1 or greater must be installed in a way that HFold can find it.    
To test if your Mac or Linux system already has CMake, you can type into a terminal:      
```
cmake --version
```
If it does not print a cmake version greater than or equal to 3.1, you will have to install CMake depending on your operating system.

#### Mac:    
Easiest way is to install homebrew and use that to install CMake.    
To do so, run the following from a terminal to install homebrew:      
```  
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"   
```    
When that finishes, run the following from a terminal to install CMake.     
```   
brew install cmake   
``` 
#### Linux:    
Run from a terminal     
```
wget http://www.cmake.org/files/v3.8/cmake-3.8.2.tar.gz
tar xzf cmake-3.8.2.tar.gz
cd cmake-3.8.2
./configure
make
make install
```
[Linux instructions source](https://geeksww.com/tutorials/operating_systems/linux/installation/downloading_compiling_and_installing_cmake_on_linux.php)

#### Steps for installation   
1. [Download the repository](https://github.com/HosnaJabbari/HFold_iterative.git) and extract the files onto your system.
2. From a command line in the root directory (where this README.md is) run
```
cmake -H. -Bbuild
cmake --build build
```   
If you need to specify a specific compiler, such as g++, you can instead run something like   
```
cmake -H. -Bbuild -DCMAKE_CXX_COMPILER=g++
cmake --build build
```   
This can be useful if you are getting errors about your compiler not having C++11 features.

Help
========================================

```
Usage: Iterative-HFold[options] [input sequence]
```

Read input file from cmdline; predict minimum free energy and optimum structure using the RNA folding algorithm.


```
  -h, --help             Print help and exit
  -V, --version          Print version and exit
  -v, --verbose          Give verbose output
  -r, --input-structure  Give a restricted structure as an input structure
  -i, --input-file       Give a path to an input file containing the sequence (and input structure if known)
  -o, --output-file      Give a path to an output file which will the sequence, and its structure and energy
  -n, --opt              Specify the number of suboptimal structures to output (default is 1)
  -d  --dangles          Specify the dangle model to be used
  -P, --paramFile        Read energy parameters from paramfile, instead of using the default parameter set.\n
      --noConv           Do not convert DNA into RNA. This will use the Matthews 2004 parameters for DNA
```

#### How to use:

        Remarks:
            make sure the <arguments> are enclosed in "", for example -r "..(...).." instead of -r ..(...)..
            The sequence does not need to be enclosed and can be given before or after the other arguments
            if no structure is provided through -r , the input structure will be the hotspot with the lowest free energy
            if -i is provided with just a file name without a path, it is assuming the file is in the diretory where the executable is called
            if -o is provided with just a file name without a path, the output file will be generated in the diretory where the executable is called
            if -v is provided, a verbose output will be given (method used is outputted)
            if -V is provided, the version is given
            if -n is provided with a number, it will modify the number of hotspots looked and outputs given (the base is 1), repeated structures are skipped. That is, if different input structures come to the same conclusion, only those that are different are shown
            If no input structure is given, or suboptimal structures are greater than the number given, CParty generates hotspots to be used as input structures -- where hotspots are energetically favorable stems
            The default parameter file is Turner2004. This can be changed via -P and specifying the parameter file you would like
    
    Sequence requirements:
        containing only characters GCAU

    Structure requirements:
        -pseudoknot free
        -containing only characters .x()
        Remarks:
            Restricted structure symbols:
                () restricted base pair
                . no restriction
                x restricted to unpaired

    Input file requirements:
            Line1: Name (optional, but must be fasta format; ignored in final input)
            Line2: Sequence (required)
            Line3: Structure (optional)
        sample:
            >Srp_005
            GCAACGAUGACAUACAUCGCUAGUCGACGC
            (............................)

#### Example:
    assume you are in the directory where the HFold_iterative executable is loacted
    ./build/Iterative-HFold -i "/home/username/Desktop/myinputfile.txt"
    ./build/Iterative-HFold -i "/home/username/Desktop/myinputfile.txt" -o "outputfile.txt"
    ./build/Iterative-HFold -i "/home/username/Desktop/myinputfile.txt" -o "/home/username/Desktop/some_folder/outputfile.txt"
    ./build/Iterative-HFold GCAACGAUGACAUACAUCGCUAGUCGACGC -r "(____________________________)"
    ./build/Iterative-HFold GCAACGAUGACAUACAUCGCUAGUCGACGC -r "(____________________________)" -o "outputfile.txt"
    ./build/Iterative-HFold -d1 GCAACGAUGACAUACAUCGCUAGUCGACGC
    ./build/Iterative-HFold GCAACGAUGACAUACAUCGCUAGUCGACGC -n 10
    ./build/Iterative-HFold -P "src/params/parameters_DP09.txt" GCAACGAUGACAUACAUCGCUAGUCGACGC 

    
## Questions
For questions, you can email mateo2@ualberta.ca
