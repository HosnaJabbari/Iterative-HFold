# make a directory in $PREFIX to build in
mkdir -p $PREFIX/Iterative-HFold && mv -v * $PREFIX/Iterative-HFold
cd $PREFIX/Iterative-HFold
# build
cmake -H. -Bbuild -DCMAKE_CXX_COMPILER=g++
cmake --build build
# make a bin folder and copy the HFold binary into it
# mkdir -p $PREFIX/bin && cp Knotty $PREFIX/bin/Knotty
# optional: removes all non-directory files 
find . -maxdepth 1 -type f -delete