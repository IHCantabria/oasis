
module purge
ml Armadillo/10.5.3-foss-2020b
ml OpenFAST/3.0.0-foss-2020b
ml OpenMPI/4.0.5-GCC-10.2.0

rm ./bin/oasis
# rm -rf build
cmake -B./build/ -H./
# cmake -DCMAKE_BUILD_TYPE=Release -B./build/ -H./
cmake --build ./build
cp ./build/OASIS ./bin/oasis