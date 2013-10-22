cd build_master
./searcher_tester \
--coarse_vocabs_file="../sift1M_double_4096.dat" \
--fine_vocabs_file="../sift1M_double_4096_8.dat" \
--query_point_type="BVEC" \
--use_residuals \
--space_dim=128 \
--subspaces_centroids_count=1024 \
--index_files_prefix="/sata/ResearchData/BigAnn/indices/sift1M_double_4096_8" \
--queries_file="/sata/ResearchData/BigAnn/bases/sift1B_queries.bvecs" \
--groundtruth_file="/sata/ResearchData/BigAnn/gnd/sift1M_groundtruth.ivecs" \
--queries_count=500 \
--neighbours_count=10000 \
--report_file="sift1M_4096_8_report.txt" \
--do_rerank
