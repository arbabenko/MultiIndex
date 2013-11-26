mkdir -p build_localVoc
cd build_localVoc
rm ./CMakeCache.txt
cmake -DCMAKE_BUILD_TYPE=Debug ..
make