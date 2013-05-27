mkdir -p build_pca_naive
cd build_pca_naive
rm ./CMakeCache.txt
cmake -DCMAKE_BUILD_TYPE=Debug ..
make