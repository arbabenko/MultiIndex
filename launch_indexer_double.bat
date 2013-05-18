cd build/Release/
indexer_launcher.exe ^
--threads_count=20 ^
--multiplicity=2 ^
--points_file="C:\Downloads\sift\sift_base.fvecs" ^
--coarse_vocabs_file="../../sift1M_double_4096.dat" ^
--fine_vocabs_file="../../sift1M_double_4096_8.dat" ^
--input_point_type="FVEC" ^
--points_count=1000000 ^
--space_dim=128 ^
--files_prefix="sift1M_double_4096_8" ^
--coarse_quantization_file="sift1M_double_4096_coarse_quantizations.bin" ^
--metainfo_file="fake.txt" ^
--use_residuals ^
--build_coarse
pause