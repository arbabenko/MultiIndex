/** @file */
// Copyright 2012 Yandex Artem Babenko
#ifndef SEARCHER_H_
#define SEARCHER_H_

#include <algorithm>
#include <map>

#include <yael/vector.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>

#include <mkl_cblas.h>

#include "data_util.h"
#include "ordered_lists_merger.h"
#include "perfomance_util.h"

extern int THREADS_COUNT;

extern Dimensions SPACE_DIMENSION;

extern enum PointType point_type;

extern int coarseM;
extern int coarseK;
extern int rerankM;
extern int rerankK;

int vocPerCoarse = rerankM / coarseM;

/**
 * \typedef This typedef is used in the first stage of search when
 * we get nearest centroids for each coarse subpace
 */
typedef vector<pair<Distance, ClusterId> > NearestSubspaceCentroids;

/**
 * This is the main class for nearest neighbour search using multiindex
 */
template<class Record, class MetaInfo>
class MultiSearcher {
 public:
 /**
  * Default constructor
  */
   MultiSearcher(int multiplicity) : multiplicity_(multiplicity);
 /**
  * Initiation function
  * @param index_files_prefix prefix of multiindex files providing the search
  * @param coarse_vocabs_filename file with coarse vocabs
  * @param fine_vocabs_filename file with fine vocabs for reranking
  * @param mode reranking approach
  * @param do_rerank should algorithm rerank short list or not
  */
  void Init(const string& index_filename,
            const string& cell_edges_filename
            const string& coarse_rotation_filename,
            const string& rerank_rotations_filename,
            const string& coarse_vocabs_filename,
            const string& rerank_vocabs_filename,
            const RerankMode& mode,
            const int subspace_centroids_to_consider,
            bool do_rerank);
 /**
  * Main interface function
  * @param point query point
  * @param k number of neighbours to get
  * @param subpace_centroids_to_consider it defines the size of working index table
  * @param neighbours result - vector of point identifiers ordered by increasing of distance to query
  */
  void GetNearestNeighbours(const Point& point, int k, 
                            vector<pair<Distance, MetaInfo> >* neighbours) const;
 /**
  * Returns searcher perfomance tester
  */
  PerfTester& GetPerfTester();
 private:
 /**
  * Function gets some nearest centroids for each coarse subspace
  * @param point query point
  * @param subspace_centroins_count how many nearest subcentroids to get
  * @param subspaces_short_lists result
  */
  void GetNearestSubspacesCentroids(const Point& point,
                                    const int subspace_centroins_count,
                                    vector<NearestSubspaceCentroids>* subspaces_short_lists) const;

 /**
  * This fuctions traverses another cell of multiindex table 
  * @param point query point
  * @param nearest_subpoints vector algorithm adds nearest neighbours in
  */
  bool TraverseNextMultiIndexCell(const Point& point,
                                  vector<pair<Distance, MetaInfo> >* nearest_subpoints) const;
 /**
  * This fuctions converts cells coordinates to appropriate range in array 
  * @param cell_coordinates coordinates of the cell
  * @param cell_start first index of range
  * @param cell_finish last index of range
  */
inline void GetCellEdgesInMultiIndexArray(const vector<int>& cell_coordinates,
                                          int* cell_start, int* cell_finish) const;
  int multiplicity_;
 /**
  * This fuctions converts complex objects to arrays and
  * pointers for usage in BLAS
  */
  void PrecomputeData();
 /**
  * Lists of coarse centroids
  */
  vector<float*> coarse_vocabs_;
 /**
  * Lists of fine centroids
  */
  vector<float*> fine_vocabs_;
 /**
  * Merger for ordered merging subspaces centroids lists
  */
  mutable OrderedListsMerger<Distance, ClusterId> merger_;
 /**
  * Should algorithm use reranking or not
  */
  bool do_rerank_;
 /**
  * Searcher perfomance tester
  */
  mutable PerfTester perf_tester_;
 /**
  * Common prefix of every index files
  */
  string index_files_prefix_;
 /**
  * Multiindex data structures
  */
  MultiIndex<Record> multiindex_;
 /**
  * Reranking approach
  */
  RerankMode rerank_mode_;
 /**
  * Struct for BLAS
  */
  vector<vector<float> > coarse_centroids_norms_;
 /**
  * Struct for BLAS
  */
  mutable Coord* products_;
 /**
  * Struct for BLAS
  */
  mutable float* residual_;
 /**
  * Number of nearest to query centroids
  * to consider for each dimension
  */
  int subspace_centroids_to_consider_;
 /**
  * Number of neighbours found to this moment
  */
  mutable int found_neghbours_count_;
 /**
  * Rotation of initial vectors (OPQ)
  */
  float* coarse_rotation_;
 /**
  * Rotation of residuals (OPQ)
  */
  vector<float*> residuals_rotations_;

};

template<class Record, class MetaInfo>
inline void RecordToMetainfoAndDistance(const Coord* point,
                                        const Record& record,
                                        pair<Distance, MetaInfo>* result,
                                        const vector<int>& cell_coordinates,
                                        const vector<float*>& fine_vocabs) {
}

/////////////// IMPLEMENTATION /////////////////////

template<class Record, class MetaInfo>
MultiSearcher<Record, MetaInfo>::MultiSearcher() {
}

template<class T>
inline T readFromFile(std::ifstream& input) {
    T result;
    input.read((char*)&result, sizeof(T));
    return result;
}

template<>
inline RerankADC8 readFromFile(std::ifstream& input) {
    RerankADC8 record;
    input.read((char*)&(record.pid), sizeof(int));
    for(q = 0; q < 8; ++q){
        input.read((char*)&(record.quantizations[q]), sizeof(unsigned char));
    }
    return record;
}

template<>
inline RerankADC16 readFromFile(std::ifstream& input) {
    RerankADC16 record;
    input.read((char*)&(record.pid), sizeof(int));
    for(q = 0; q < 16; ++q){
        input.read((char*)&(record.quantizations[q]), sizeof(unsigned char));
    }
    return record;
}


template<class T>
void readSubVector(const string& filename,
                   int start_record_id,
                   int chunk_size,
                   int records_count,
                   vector<T>* output) {
    chunk_size = std::min(records_count - start_record_id, chunk_size);
    std::ifstream input(filename.c_str(), ios::binary | ios::in);
    input.seekg(start_record_id * sizeof(T));
    for(int i = 0; i < chunk_size; ++i) {
      output->at(start_record_id + i) = readFromFile(input);
    }
}

template<class T>
void fillVector<T>(const string& filename,
                   vector<T>* output) {
  std::ifstream in(filename, std::ifstream::in | std::ifstream::binary);
  in.seekg(0, std::ifstream::end);
  records_count = in.tellg() / sizeof(T);
  in.close();
  output->resize(records_count);
  int chunk_size = records_count / THREADS_COUNT;
  boost::thread_group index_threads;
  for(int thread_id = 0; thread_id < THREADS_COUNT; ++thread_id) {
    int start_record_id = chunk_size * thread_id;
    index_threads.create_thread(boost::bind(&readSubVector, filename, start_record_id,
                                            chunk_size, output));
  }
  index_threads.join_all();
}

template<class Record, class MetaInfo>
void MultiSearcher<Record, MetaInfo>::Init(const string& index_filename,
                                           const string& cell_edges_filename
                                           const string& coarse_rotation_filename,
                                           const string& rerank_rotations_filename,
                                           const string& coarse_vocabs_filename,
                                           const string& rerank_vocabs_filename,
                                           const RerankMode& mode,
                                           const int subspace_centroids_to_consider,
                                           bool do_rerank) {
  do_rerank_ = do_rerank;
  subspace_centroids_to_consider_ = subspace_centroids_to_consider;
  rerank_mode_ = mode;
  merger_.GetYieldedItems().table.resize(std::pow((float)subspace_centroids_to_consider,
		                                         (int)coarse_vocabs_.size()));
  for(int i = 0; i < multiplicity_; ++i) {
    merger_.GetYieldedItems().dimensions.push_back(subspace_centroids_to_consider);
  }
  fvecs_fread(fopen(coarse_rotation_filename, "r"), coarse_rotation_, dim, dim);
  FILE* rerank_matrix_file = fopen(rerank_rotations_filename, "r");
  residuals_rotations_.resize(multiplicity_);
  for (int m = 0; m < multiplicity_; ++m){
    fvecs_fread(rerank_matrix_file, residuals_rotations_[m], 
                SPACE_DIMENSION / multiplicity_, SPACE_DIMENSION / multiplicity_);
  }
  fclose(rerank_matrix_file);
  FILE* coarse_vocabs_file = fopen(coarse_vocabs_filename, "r");
  coarse_vocabs_.resize(multiplicity_);
  for (int m = 0; m < multiplicity_; ++m){
    fvecs_fread(coarse_vocabs_file, coarse_vocabs_[m], coarseK, SPACE_DIMENSION / multiplicity_);
  }
  fclose(coarse_vocabs_file);
  FILE* rerank_vocabs_file = fopen(rerank_vocabs_filename, "r");
  fine_vocabs_.resize(rerankM);
  for (int m = 0; m < rerankM; ++m){
    fvecs_fread(rerank_vocabs_file, fine_vocabs_[m], rerankK, SPACE_DIMENSION / rerankM);
  }
  fclose(rerank_vocabs_file);
  multiindex_.cell_edges.length = std::pow((float)coarseK, coarseM);
  fillVector<int>(cell_edges_filename, &(multiindex_.cell_edges.table));
  fillVector<Record>(index_filename, &(multiindex_.multiindex));
  PrecomputeData();
}

template<class Record, class MetaInfo>
void MultiSearcher<Record, MetaInfo>::PrecomputeData(){
  coarse_centroids_norms_.resize(coarseM, coarseK);
  int vocab_dim = SPACE_DIMENSION / coarseM;
  for(int coarse_id = 0; coarse_id < coarseM; ++coarse_id) {
    for(int k = 0; k < coarseK; ++k) {
      coarse_centroids_norms_[coarse_id][k] = cblas_sdot(vocab_dim, 
                                                         coarse_vocabs_[coarse_id] + vocab_dim * k, 1,
                                                         coarse_vocabs_[coarse_id] + vocab_dim * k, 1);
    }
  }
  products_ = new Coord[coarseK];
  residual_ = new Coord[dim];
}

template<class Record, class MetaInfo>
PerfTester& MultiSearcher<Record, MetaInfo>::GetPerfTester() {
  return perf_tester_;
}

template<class Record, class MetaInfo>
void MultiSearcher<Record, MetaInfo>::GetNearestSubspacesCentroids(const Point& point,
                                                                   const int subspace_centroins_count,
                                                                   vector<NearestSubspaceCentroids>*
                                                                   subspaces_short_lists) const {
  subspaces_short_lists->resize(coarse_vocabs_.size());
  Dimensions subspace_dimension = point.size() / coarse_vocabs_.size();
  for(int subspace_index = 0; subspace_index < coarse_vocabs_.size(); ++subspace_index) {
    Dimensions start_dim = subspace_index * subspace_dimension;
    Dimensions final_dim = std::min((Dimensions)point.size(), start_dim + subspace_dimension);
    memset(products_, 0, sizeof(float) * coarseK);
    cblas_saxpy(coarseK, 1, &(coarse_centroids_norms_[subspace_index][0]), 1, products_, 1);
    cblas_sgemv(CblasRowMajor, CblasNoTrans, coarseK, subspace_dimension, -2.0,
                coarse_vocabs_[subspace_index], subspace_dimension, &(point[start_dim]), 1, 1, products_, 1);
    subspaces_short_lists->at(subspace_index).resize(coarseK);
    for(int i = 0; i < coarseK; ++i) {
      subspaces_short_lists->at(subspace_index)[i] = std::make_pair(products_[i], i);
    }
    std::nth_element(subspaces_short_lists->at(subspace_index).begin(),
                     subspaces_short_lists->at(subspace_index).begin() + subspace_centroins_count,
                     subspaces_short_lists->at(subspace_index).end());
    subspaces_short_lists->at(subspace_index).resize(subspace_centroins_count);
    std::sort(subspaces_short_lists->at(subspace_index).begin(),
              subspaces_short_lists->at(subspace_index).end());
  }
}

template<class Record, class MetaInfo>
void MultiSearcher<Record, MetaInfo>::GetCellEdgesInMultiIndexArray(const vector<int>& cell_coordinates,
                                                                    int* cell_start, int* cell_finish) const {
  int global_index = multiindex_.cell_edges.GetCellGlobalIndex(cell_coordinates);
  *cell_start = multiindex_.cell_edges.table[global_index];
  if(global_index + 1 == multiindex_.cell_edges.table.size()) {
    *cell_finish = multiindex_.multiindex.size();
  } else {
    *cell_finish = multiindex_.cell_edges.table[global_index + 1];
  }
}

template<class Record, class MetaInfo>
bool MultiSearcher<Record, MetaInfo>::TraverseNextMultiIndexCell(const Point& point,
                                                                 vector<pair<Distance, MetaInfo> >*
                                                                             nearest_subpoints) const {
  MergedItemIndices cell_inner_indices;
  clock_t before = clock();
  if(!merger_.GetNextMergedItemIndices(&cell_inner_indices)) {
    return false;
  }
  clock_t after = clock();
  perf_tester_.cell_coordinates_time += after - before;
  vector<int> cell_coordinates(cell_inner_indices.size());
  for(int list_index = 0; list_index < merger_.lists_ptr->size(); ++list_index) {
    cell_coordinates[list_index] = merger_.lists_ptr->at(list_index)[cell_inner_indices[list_index]].second;
  }
  int cell_start, cell_finish;
  before = clock();
  GetCellEdgesInMultiIndexArray(cell_coordinates, &cell_start, &cell_finish);
  after = clock();
  perf_tester_.cell_edges_time += after - before;
  if(cell_start >= cell_finish) {
    return true;
  }
  typename vector<Record>::const_iterator it = multiindex_.multiindex.begin() + cell_start;
  GetResidual(point, cell_coordinates, coarse_vocabs_, residual_);
  for(int m = 0; m < coarseM; ++m) {
    cblas_sgemv(CblasRowMajor, CblasNoTrans, SPACE_DIMENSION / coarseM, SPACE_DIMENSION / coarseM, 1,
                residuals_rotations_[m], SPACE_DIMENSION / coarseM,
                &(residual_[m * SPACE_DIMENSION / coarseM]), 1, 0, &(residual_[m * SPACE_DIMENSION / coarseM]), 1);      
  }
  cell_finish = std::min((int)cell_finish, cell_start + (int)nearest_subpoints->size() - found_neghbours_count_);
  for(int array_index = cell_start; array_index < cell_finish; ++array_index) {
    if(rerank_mode_ == USE_RESIDUALS) {
      RecordToMetainfoAndDistance<Record, MetaInfo>(residual_, *it,
                                                    &(nearest_subpoints->at(found_neghbours_count_)),
                                                    cell_coordinates, fine_vocabs_);
    } else if(rerank_mode_ == USE_INIT_POINTS) {
      RecordToMetainfoAndDistance<Record, MetaInfo>(&(point[0]), *it,
                                                    &(nearest_subpoints->at(found_neghbours_count_)),
                                                    cell_coordinates, fine_vocabs_);
    }
    perf_tester_.NextNeighbour();
    ++found_neghbours_count_;
    ++it;
  }
  return true;
}


template<class Record, class MetaInfo>
void MultiSearcher<Record, MetaInfo>::GetNearestNeighbours(const Point& point, int k, 
                                                           vector<pair<Distance, MetaInfo> >* neighbours) const {
  assert(k > 0);
  perf_tester_.handled_queries_count += 1;
  neighbours->resize(k);
  perf_tester_.ResetQuerywiseStatistic();
  clock_t start = clock();
  perf_tester_.search_start = start;
  clock_t before = clock();
  // rotate point (OPQ)
  cblas_sgemv(CblasRowMajor, CblasNoTrans, SPACE_DIMENSION, SPACE_DIMENSION, 1,
              coarse_rotation_, SPACE_DIMENSION, &(point[0]), 1, 0, &(point[0]), 1);
  vector<NearestSubspaceCentroids> subspaces_short_lists;
  assert(subspace_centroids_to_consider_ > 0);
  GetNearestSubspacesCentroids(point, subspace_centroids_to_consider_, &subspaces_short_lists);
  clock_t after = clock();
  perf_tester_.nearest_subcentroids_time += after - before;
  clock_t before_merger = clock();
  merger_.setLists(subspaces_short_lists);
  clock_t after_merger = clock();
  perf_tester_.merger_init_time += after_merger - before_merger;
  clock_t before_traversal = clock();
  found_neghbours_count_ = 0;
  bool traverse_next_cell = true;
  int cells_visited = 0;
  while(found_neghbours_count_ < k && traverse_next_cell) {
    perf_tester_.cells_traversed += 1;
    traverse_next_cell = TraverseNextMultiIndexCell(point, neighbours);
    cells_visited += 1;
  }
  clock_t after_traversal = clock();
  perf_tester_.full_traversal_time += after_traversal - before_traversal;
  if(do_rerank_) {
    std::sort(neighbours->begin(), neighbours->end());
  }
  clock_t finish = clock();
  perf_tester_.full_search_time += finish - start;
}

template<>
inline void RecordToMetainfoAndDistance<RerankADC8, PointId>(const Coord* point, const RerankADC8& record,
                                                             pair<Distance, PointId>* result,
                                                             const vector<int>& cell_coordinates,
                                                             const vector<float*>& fine_vocabs) {
  result->second = record.pid;
  int subvectors_dim = SPACE_DIMENSION / rerankM;
  char* rerank_info_ptr = (char*)&record + sizeof(record.pid);
  for(int centroid_index = 0; centroid_index < rerankM; ++centroid_index) {
    int start_dim = centroid_index * subvectors_dim;
    int final_dim = start_dim + subvectors_dim;
    FineClusterId pid_nearest_centroid = *((FineClusterId*)rerank_info_ptr);
    rerank_info_ptr += sizeof(FineClusterId);
    Distance subvector_distance = 0;
    for(int i = start_dim; i < final_dim; ++i) {
      Coord diff = fine_vocabs[centroid_index][subvectors_dim * pid_nearest_centroid + i - start_dim] - point[i];
      subvector_distance += diff * diff;
    }
    result->first += subvector_distance;
  }
}

template<>
inline void RecordToMetainfoAndDistance<RerankADC16, PointId>(const Coord* point, const RerankADC16& record,
                                                              pair<Distance, PointId>* result,
                                                              const vector<int>& cell_coordinates,
                                                              const vector<float*>& fine_vocabs) {
  result->second = record.pid;
  int subvectors_dim = SPACE_DIMENSION / rerankM;
  char* rerank_info_ptr = (char*)&record + sizeof(record.pid);
  for(int centroid_index = 0; centroid_index < rerankM; ++centroid_index) {
    int start_dim = centroid_index * subvectors_dim;
    int final_dim = start_dim + subvectors_dim;
    FineClusterId pid_nearest_centroid = *((FineClusterId*)rerank_info_ptr);
    rerank_info_ptr += sizeof(FineClusterId);
    Distance subvector_distance = 0;
    for(int i = start_dim; i < final_dim; ++i) {
      Coord diff = fine_vocabs[centroid_index][subvectors_dim * pid_nearest_centroid + i - start_dim] - point[i];
      subvector_distance += diff * diff;
    }
    result->first += subvector_distance;
  }
}

template class MultiSearcher<RerankADC8, PointId>;
template class MultiSearcher<RerankADC16, PointId>;
template class MultiSearcher<PointId, PointId>;

#endif

