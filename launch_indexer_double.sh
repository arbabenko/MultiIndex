cd build_master
./indexer_launcher \
--threads_count=32 \
--multiplicity=2 \
--points_file="/sata/ResearchData/BigAnn/bases/sift1M.bvecs" \
--coarse_vocabs_file="../sift1M_double_4096.dat" \
--fine_vocabs_file="../sift1M_double_4096_8.dat" \
--input_point_type="BVEC" \
--points_count=1000000 \
--space_dim=128 \
--files_prefix="/sata/ResearchData/BigAnn/indices/sift1M_double_4096_8" \
--coarse_quantization_file="/sata/ResearchData/BigAnn/cq/sift1M_double_4096_coarse_quantizations.bin" \
--metainfo_file="fake.txt" \
--use_residuals \
--build_coarse

