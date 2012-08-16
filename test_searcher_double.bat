cd build/Release/
searcher_tester.exe ^
--coarse_vocabs_file="../../sift1M_double_4096.dat" ^
--fine_vocabs_file="../../sift1M_double_4096_8.dat" ^
--query_point_type="FVEC" ^
--use_residuals ^
--space_dim=128 ^
--subspaces_centroids_count=1024 ^
--index_files_prefix="sift1M_double_4096_8" ^
--queries_file="C:\Downloads\sift\sift_query.fvecs" ^
--groundtruth_file="C:\Downloads\sift/sift_groundtruth.ivecs" ^
--queries_count=10000 ^
--neighbours_count=10000 ^
--report_file="sift1B_16384_8_report.txt" ^
--do_rerank
pause