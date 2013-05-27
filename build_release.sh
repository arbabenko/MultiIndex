mkdir -p build_pca_smart
cd build_pca_smart
rm ./CMakeCache.txt
cmake -DCMAKE_BUILD_TYPE=Release ..
make