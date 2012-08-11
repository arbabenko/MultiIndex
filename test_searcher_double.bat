cd build/Release/
searcher_tester.exe ^
--coarse_vocabs_file="../../sift1M_double_4096.dat" ^
--fine_vocabs_file="../../sift1M_double_4096_8.dat" ^
--query_point_type="BVEC" ^
--use_residuals ^
--space_dim=128 ^
--subspaces_centroids_count=256 ^
--index_files_prefix="sift1M_double_16384_8" ^
--queries_file="../../queries.bvecs" ^
--groundtruth_file="../../groundtruth.ivecs" ^
--queries_count=100 ^
--neighbours_count=10000 ^
--report_file="sift1B_16384_8_report.txt" ^
--do_rerank
pause