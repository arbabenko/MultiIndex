#!/usr/bin/python
import os
import datetime
import sys

############# launch configuration #################

# folder to keep built binaries in
build_folder = './build_hier'

# number of threads to use (max = 32)
threads_count = 32

# Folder with BigAnn base
bigann_root = '/sata/ResearchData/BigAnn'

# input point type (BVEC or FVEC)
input_type = 'BVEC'

# file with input point (.bvecs or .fvecs)
points_file = 'sift1M.bvecs'

# prefix of all vocabs, coarse quantizations, etc.
prefix = 'sift1M'

# input points count
points_count = 1000000

# dimension of input space
space_dim = 128

# main vocabs size
main_vocabs_size = 16384

# residual vocabs count
res_vocabs_count = 16384

# should we use residuals for reranking
use_residuals = 1

# should we calculate coarse quantizations?
build_coarse = 0

# postfix added by users to all multiindex files
user_added_postfix = ''

##################################################

multiplicity_extension = 'hier'
#if multiplicity == 1:
#    multiplicity_extension = 'single'
#if multiplicity == 2:
#    multiplicity_extension = 'double'
#if multiplicity == 4:
#    multiplicity_extension = 'quad'

main_vocabs_filename = prefix + '_' + multiplicity_extension + '_' + str(main_vocabs_size) + '.dat'
res_vocabs_filename = prefix + '_' + multiplicity_extension + '_' + str(main_vocabs_size) + '_' + str(res_vocabs_count) + '.dat'
filename_prefix = prefix + '_' + multiplicity_extension + '_' + str(main_vocabs_size) + '_' + str(res_vocabs_count) + user_added_postfix
coarse_quantization_filename = prefix + '_' + multiplicity_extension + '_' + str(main_vocabs_size) + user_added_postfix + '_coarse_quantizations.bin'

launch_time = datetime.datetime.now().strftime("%I_%M%p_%B_%d_%Y")
os.system('mkdir -p ' + build_folder + '/' + launch_time)
os.system('cp ' + build_folder + '/indexer_launcher ' + build_folder + '/' + launch_time)
os.system('cp run_indexer.py ' + build_folder + '/' + launch_time)

launch_line = build_folder + '/' + launch_time + '/indexer_launcher '
launch_line = launch_line + '--threads_count=' + str(threads_count) + ' '
launch_line = launch_line + '--points_file=' + bigann_root + '/bases/' + points_file + ' '
launch_line = launch_line + '--main_vocabs_file=' + bigann_root + '/coarse_vocabs/' + main_vocabs_filename + ' '
launch_line = launch_line + '--res_vocabs_file=' + bigann_root + '/fine_vocabs/' + res_vocabs_filename + ' '
launch_line = launch_line + '--input_point_type=' + input_type + ' '
launch_line = launch_line + '--points_count=' + str(points_count) + ' '
launch_line = launch_line + '--space_dim=' + str(space_dim) + ' '
launch_line = launch_line + '--files_prefix=' + bigann_root + '/indices/' + filename_prefix + ' ' 
launch_line = launch_line + '--coarse_quantization_file=' + bigann_root + '/cq/' + coarse_quantization_filename + ' '
launch_line = launch_line + '--metainfo_file=fake.txt' + ' '
if use_residuals:
    launch_line = launch_line + '--use_residuals' + ' '
if build_coarse:
    launch_line = launch_line + '--build_coarse' + ' '

f = open(build_folder + '/' + launch_time + '/launch.sh', 'w')
f.write(launch_line)
f.close()
os.system('nohup ' + launch_line + ' > ' + build_folder + '/' + launch_time + '/log.txt' + ' &')
print 'Log file: ' + build_folder + '/' + launch_time + '/log.txt'
