<!-- Edited Aug 5 2022 by Connie He to reflect conda package installation-->
# Iterative HFold

#### Description:
Software implementation of Iterative HFold.      
Iterative HFold is an algorithm for predicting the pseudoknotted secondary structures of RNA using relaxed Hierarchical Folding. 

Paper: https://www.researchgate.net/publication/262810273_A_fast_and_robust_iterative_algorithm_for_prediction_of_RNA_pseudoknotted_secondary_structures

On the dataset tested in this paper, Iterative HFold generally has better accuracy that its predecessor, [HFold](https://github.com/HosnaJabbari/HFold).

#### Supported OS: 
Linux, macOS

### Installation:
Requirements: conda.    

You can get conda by installing either [Miniconda](https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links) or [Anaconda](https://www.anaconda.com/products/distribution), but the easier method is to install Miniconda.

To test if your Mac or Linux system already has conda, you can type into a terminal:      
```
conda --version
```
If it prints the message, `conda: command not found`, you will have to install conda. 

#### Conda Installation:
In the terminal, run one of the following set of commands according to your operating system:

##### Mac (Intel chip):
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

##### Mac (M1 chip):
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
bash Miniconda3-latest-MacOSX-arm64.sh
```

##### Linux 64-bit:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

##### Linux (others)
For other Linux platforms, use the previous Linux 64-bit commands but replace `x86_64` with one of the following as per your platform:

- AWS Graviton2 / ARM64: `aarch64`

- Power8 / Power9: `ppc64le`

- IBM Z / LinuxOne: `s390x`

<br>

After running one of the previous sets of commands, follow the prompts from the installer. When prompted with 
```
Do you wish the installer to initialize Miniconda3
by running conda init? [yes|no]
```
Enter `yes`. 
Once the installation is complete, restart your terminal and run the following in the terminal to verify that conda was installed:     
```   
conda --version  
``` 
If it does not work, run `conda init --all` then restart your terminal.

#### Steps for installation   
To install Iterative-HFold, run the following command in the terminal:
```
conda install -c uvic-cobra iterative-hfold
```
And enter `y` when prompted with `Proceed ([y]/n)?`.

#### How to use:
    Arguments:
        Iterative-HFold:
            --s <sequence>
            --r <structure>
            --i </path/to/file>
            --o </path/to/file>

        Remarks:
            make sure the <arguments> are enclosed in "", for example --r "..().." instead of --r ..()..
            input file for --i must be .txt
            if --i is provided with just a file name without a path, it is assuming the file is in the diretory where the executable is called
            if --o is provided with just a file name without a path, the output file will be generated in the diretory where the executable is called
            if --o is provided with just a file name without a path, and if --i is provided, then the output file will be generated in the directory where the input file is located
    
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

#### Example:
    Iterative-HFold --i "/home/username/Desktop/myinputfile.txt"
    Iterative-HFold --i "/home/username/Desktop/myinputfile.txt" -o "outputfile.txt"
    Iterative-HFold --i "/home/username/Desktop/myinputfile.txt" -o "/home/username/Desktop/some_folder/outputfile.txt"
    Iterative-HFold --s "GCAACGAUGACAUACAUCGCUAGUCGACGC" -r "(____________________________)"
    Iterative-HFold --s "GCAACGAUGACAUACAUCGCUAGUCGACGC" -r "(____________________________)" -o "outputfile.txt"
    Iterative-HFold --s "GCAACGAUGACAUACAUCGCUAGUCGACGC" -r "(____________________________)" -o "/home/username/Desktop/some_folder/outputfile.txt"

    
#### Exit code:
    0       success
    1	    invalid argument error 
    3	    thread error
    4       i/o error
    5       pipe error
    6       positive energy error
    error code with special meaning: http://tldp.org/LDP/abs/html/exitcodes.html
    2	    Misuse of shell builtins (according to Bash documentation)
    126	    Command invoked cannot execute
    127	    "command not found"
    128	    Invalid argument to exit	
    128+n	Fatal error signal "n"
    130	    Script terminated by Control-C
    255	    Exit status out of range (range is 0-255)