(Works on linux and macOS)

## Install build tools
- Download conda and run `conda install anaconda-client conda-build` to install the build tools and upload tools

## How to build the package
- Start in the Iterative-HFold  folder
- Run `conda build ./conda_recipe` to build package
- Run `conda install --use-local Iterative-HFold` to install (case sensitive)
- Run `Iterative-HFold  --help` to insure it installed properly (case sensitive)

## Setup Conda Account (For uploading packages)
- Create an account at https://anaconda.org/
- Type `anaconda login` in terminal and login to your account
- Run `conda config --set anaconda_upload no` to prevent automatic uploads to anaconda after every build

## Upload package
- Find the package you built using `conda build ./conda_recipe --output`
- Type `anaconda upload -u COBRALab /path/to/package.conda`
- You must be added to the organization for this to work
- Please update the version in the meta.yaml file before uploading