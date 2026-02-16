module purge
ml Armadillo/10.5.3-foss-2020b
ml FASTurbine_wrapper/0.1.1-foss-2020b
ml stl_reader/1.0.0-foss-2020b

rm -rf bin || true
mkdir bin

# rm -rf build
# cmake -DCMAKE_BUILD_TYPE=Debug -DOASIS_USE_OPENFAST=ON -B./build/ -H./
cmake -DCMAKE_BUILD_TYPE=Release -DOASIS_USE_OPENFAST=ON -B./build/ -H./

cmake --build ./build

cp ./build/OASIS ./bin/oasis
