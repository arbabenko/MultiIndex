cd build_localVoc
./searcher_tester \
--index_file="/sata/ResearchData/indexLocal_1000000000.dat" \
--cell_edges="/sata/ResearchData/cellStarts_1000000000.dat" \
--coarse_rotation_file="/sata/ResearchData/coarseRotation.fvecs" \
--rerank_rotation_file="/sata/ResearchData/rerankRotations.fvecs" \
--coarse_vocabs_file="/sata/ResearchData/coarseVoc.fvecs" \
--rerank_vocabs_file="/sata/ResearchData/rerankVocabsLocal.fvecs" \
--query_point_type="BVEC" \
--use_residuals \
--space_dim=128 \
--subspaces_centroids_count=1024 \
--queries_file="/sata/ResearchData/BigAnn/bases/sift1B_queries.bvecs" \
--groundtruth_file="/sata/ResearchData/BigAnn/gnd/sift1B_groundtruth.ivecs" \
--queries_count=1000 \
--neighbours_count=10000 \
--report_file="sift1M_4096_8_report.txt" \
--multi=2 \
--do_rerank \
--local_vocabs_mode 




