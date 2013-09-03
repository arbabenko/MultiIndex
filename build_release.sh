mkdir -p build_hier
cd build_hier
rm ./CMakeCache.txt
cmake -DCMAKE_BUILD_TYPE=Release ..
make
