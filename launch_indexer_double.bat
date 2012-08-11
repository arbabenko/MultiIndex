cd build/Release/
indexer_launcher.exe ^
--threads_count=8 ^
--multiplicity=2 ^
--points_file="../../sift1M.bvecs" ^
--coarse_vocabs_file="../../sift1M_double_4096.dat" ^
--fine_vocabs_file="../../sift1M_double_4096_8.dat" ^
--input_point_type="BVEC" ^
--points_count=1000000 ^
--space_dim=128 ^
--files_prefix="sift1M_double_16384_8" ^
--coarse_quantization_file="sift1M_double_16384_coarse_quantizations.bin" ^
--metainfo_file="fake.txt" ^
--use_residuals ^
--build_coarse
pause