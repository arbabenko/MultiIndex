/** @file */
// Copyright 2012 Yandex Artem Babenko
#ifndef SEARCHER_H_
#define SEARCHER_H_

#include <algorithm>
#include <map>

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

/**
 * This is the main class for nearest neighbour search using hierarchical index
 */
template<class Record, class MetaInfo>
class Searcher {
 public:
 /**
  * Default constructor
  */
  Searcher();
 /**
  * Initiation function
  * @param index_files_prefix prefix of multiindex files providing the search
  * @param main_vocabs_filename file with vocabs for initial points
  * @param res_vocabs_filename file with vocabs for residuals
  * @param mode reranking approach
  * @param do_rerank should algorithm rerank short list or not
  */
  void Init(const string& index_files_prefix,
            const string& main_vocabs_filename,
            const string& res_vocabs_filename,
            const RerankMode& mode,
            const int main_centroids_to_consider,
            bool do_rerank);
 /**
  * Main interface function
  * @param point query point
  * @param k number of neighbours to get
  * @param main_centroids_to_consider number of main clusters to search in
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
  * This functions deserializes all structures for search
  * @param index_files_prefix prefix of multiindex files providing the search
  * @param main_vocabs_filename file with centroids for initial points
  * @param res_vocabs_filename file with centroids for residuals
  */
  void DeserializeData(const string& index_files_prefix,
                       const string& main_vocabs_filename,
                       const string& res_vocabs_filename);
 /**
  * Function gets some nearest centroids for each coarse subspace
  * @param point query point
  * @param subspace_centroins_count how many nearest subcentroids to get
  * @param subspaces_short_lists result
  */
  //void GetNearestSubspacesCentroids(const Point& point,
  //                                  const int subspace_centroins_count,
  //                                  vector<NearestSubspaceCentroids>* subspaces_short_lists) const;

 /**
  * This fuctions traverses another cell of multiindex table 
  * @param point query point
  * @param nearest_subpoints vector algorithm adds nearest neighbours in
  */
  //bool TraverseNextMultiIndexCell(const Point& point,
  //                                vector<pair<Distance, MetaInfo> >* nearest_subpoints) const;
 /**
  * This fuctions converts cells coordinates to appropriate range in array 
  * @param cell_coordinates coordinates of the cell
  * @param cell_start first index of range
  * @param cell_finish last index of range
  */
  inline void GetCellEdgesInMultiIndexArray(const vector<int>& cell_coordinates,
                                            int* cell_start, int* cell_finish) const;
 /**
  * This fuctions converts complex objects to arrays and
  * pointers for usage in BLAS
  */
  void InitBlasStructures();
 /**
  * Lists of coarse centroids
  */
  Centroids main_vocabs_;
 /**
  * Lists of fine centroids
  */
  Centroids res_vocabs_;
 /**
  * Merger for ordered merging subspaces centroids lists
  */
  //mutable OrderedListsMerger<Distance, ClusterId> merger_;
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
  //vector<float*> coarse_vocabs_matrices_;
 /**
  * Struct for BLAS
  */
  //vector<vector<float> > coarse_centroids_norms_;
 /**
  * Struct for BLAS
  */
  //mutable Coord* products_;
 /**
  * Struct for BLAS
  */
  //mutable vector<Coord> query_norms_;
 /**
  * Struct for BLAS
  */
  //mutable float* residual_;
 /**
  * Number of nearest to query centroids
  * to consider for each dimension
  */
  int main_centroids_to_consider_;
 /**
  * Number of neighbours found to this moment
  */
  mutable int found_neghbours_count_;

  vector<vector<float> > precomputed_norms_;

  vector<float> main_norms_;

  float* main_vocabs_matrix_;

  float* res_vocabs_matrix_;
};

template<class Record, class MetaInfo>
inline void RecordToMetainfoAndDistance(const Coord* point,
                                        const Record& record,
                                        pair<Distance, MetaInfo>* result,
                                        const vector<int>& cell_coordinates,
                                        const vector<Centroids>& fine_vocabs) {
}

/////////////// IMPLEMENTATION /////////////////////

template<class Record, class MetaInfo>
Searcher<Record, MetaInfo>::Searcher() {
}

template<class Record, class MetaInfo>
void Searcher<Record, MetaInfo>::DeserializeData(const string& index_files_prefix,
                                                 const string& main_vocabs_file,
                                                 const string& res_vocabs_file) {
  cout << "Data deserializing started...\n";
  ifstream cell_edges(string(index_files_prefix + "_cell_edges.bin").c_str(), ios::binary);
  if(!cell_edges.good()) {
    throw std::logic_error("Bad input cell edges stream");
  }
  boost::archive::binary_iarchive arc_cell_edges(cell_edges);
  arc_cell_edges >> multiindex_.cell_edges;
  cout << "Cell edges deserialized...\n";

  ifstream multi_array(string(index_files_prefix + "_multi_array.bin").c_str(), ios::binary);
  if(!multi_array.good()) {
    throw std::logic_error("Bad input multiarray stream");
  }
  boost::archive::binary_iarchive arc_multi_array(multi_array);
  arc_multi_array >> multiindex_.multiindex;
  cout << "Multiarray deserialized...\n";

  ifstream norms(string(index_files_prefix + "_prec_norms.bin").c_str(), ios::binary);
  if(!norms.good()) {
    throw std::logic_error("Bad input norms stream");
  }
  boost::archive::binary_iarchive arc_norms(norms);
  arc_norms >> precomputed_norms_;
  cout << "Effective norms deserialized...\n";

  ifstream main_norms(string(index_files_prefix + "_main_norms.bin").c_str(), ios::binary);
  if(!main_norms.good()) {
    throw std::logic_error("Bad input main norms stream");
  }
  boost::archive::binary_iarchive arc_main_norms(main_norms);
  arc_main_norms >> main_norms_;  
  cout << "Main norms deserialized...\n";

  ReadCentroids<float>(main_vocabs_file, SPACE_DIMENSION, &main_vocabs_);
  ReadCentroids<float>(res_vocabs_file, SPACE_DIMENSION, &res_vocabs_);

}

template<class Record, class MetaInfo>
void Searcher<Record, MetaInfo>::Init(const string& index_files_prefix,
                                      const string& main_vocabs_filename,
                                      const string& res_vocabs_filename,
                                      const RerankMode& mode,
                                      const int main_centroids_to_consider,
                                      const bool do_rerank) {
  do_rerank_ = do_rerank;
  index_files_prefix_ = index_files_prefix;
  main_centroids_to_consider_ = main_centroids_to_consider;
  DeserializeData(index_files_prefix, main_vocabs_filename, res_vocabs_filename);
  rerank_mode_ = mode;
  InitBlasStructures();
}

template<class Record, class MetaInfo>
void Searcher<Record, MetaInfo>::InitBlasStructures(){

  main_vocabs_matrix_ = new float[main_vocabs_.size() * main_vocabs_[0].size()];
  for(int cid = 0; cid < main_vocabs_.size(); ++cid) {
    for(int coord = 0; coord < main_vocabs_[0].size(); ++coord) {
        main_vocabs_matrix_[cid * main_vocabs_[0].size() + coord] = main_vocabs_[cid][coord];
    }
  }
  res_vocabs_matrix_ = new float[res_vocabs_.size() * res_vocabs_[0].size()];
  for(int cid = 0; cid < res_vocabs_.size(); ++cid) {
    for(int coord = 0; coord < res_vocabs_[0].size(); ++coord) {
        res_vocabs_matrix_[cid * res_vocabs_[0].size() + coord] = res_vocabs_[cid][coord];
    }
  }
}

template<class Record, class MetaInfo>
PerfTester& Searcher<Record, MetaInfo>::GetPerfTester() {
  return perf_tester_;
}

template<class Record, class MetaInfo>
void Searcher<Record, MetaInfo>::GetCellEdgesInMultiIndexArray(const vector<int>& cell_coordinates,
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
void Searcher<Record, MetaInfo>::GetNearestNeighbours(const Point& point, int k, 
                                                           vector<pair<Distance, MetaInfo> >* neighbours) const {
  found_neghbours_count_ = 0;
  neighbours->resize(k);
  perf_tester_.ResetQuerywiseStatistic();
  clock_t start = clock();
  perf_tester_.search_start = start;
  clock_t before = clock();
  vector<float> main_products(main_vocabs.size(), 0);
  vector<float> res_products(res_vocabs.size(), 0);
  vector<pair<float, int> > distance_to_clusterId(main_products.size());
  cblas_sgemv(CblasRowMajor, CblasNoTrans, main_vocabs.size(), main_vocabs[0].size(), 1.0,
              main_vocabs_matrix_, main_vocabs[0].size(), &(point[0]), 1, 1, &(main_products[0]), 1);
  for (int main_cid = 0; main_cid < main_vocabs.size(); ++main_cid) {
    distance_to_clusterId[main_cid] = std::make_pair(main_norms_[main_cid] / 2 - main_products[main_cid], main_cid);
  }
  cblas_sgemv(CblasRowMajor, CblasNoTrans, res_vocabs.size(), res_vocabs[0].size(), 1.0,
              res_vocabs_matrix_, res_vocabs[0].size(), &(point[0]), 1, 1, &(res_products[0]), 1);
  std::sort(distance_to_clusterId.begin(), distance_to_clusterId.end());
  clock_t after = clock();
  perf_tester_.nearest_subcentroids_time += after - before;
  clock_t before_merger = clock();
  std::multimap<float, pair<ClusterId, ClusterId> > queue;
  for(int main_cid = 0; main_cid < main_centroids_to_consider_; ++main_cid) {
    for(int res_cid = 0; res_cid < res_vocabs_.size(); ++res_cid) {
      int main_cluster = distance_to_clusterId[main_cid].second;
      float score = precomputed_norms_[main_cluster][res_cid] -
                    2 * main_products[main_cluster] -
                    2 * res_products[res_cid];
      queue.insert(std::make_pair(score, std::make_pair(main_cluster,res_cid)));
    }
  }
  clock_t after_merger = clock();
  perf_tester_.merger_init_time += after_merger - before_merger;
  cout << "Traverse started" << endl;
  std::multimap<float, pair<ClusterId, ClusterId> >::iterator current_cell = queue.begin();
  clock_t before_traversal = clock();
  while(found_neghbours_count_ < k && current_cell != queue.end()) {
      vector<ClusterId> quantization;
      quantization.push_back(current_cell->second.first);
      quantization.push_back(current_cell->second.second);
      int cell_start, cell_finish;
      GetCellEdgesInMultiIndexArray(quantization, &cell_start, &cell_finish);
      if(cell_start >= cell_finish) {
         ++current_cell;
         continue;
      }
      typename vector<Record>::const_iterator it = multiindex_.multiindex.begin() + cell_start;
      cell_finish = std::min((int)cell_finish, cell_start + (int)neighbours->size() - found_neghbours_count_);
      for(int array_index = cell_start; array_index < cell_finish; ++array_index) {
        Record record = multiindex_.multiindex[array_index];
        neighbours->at(found_neghbours_count_).second = record.pid;
        ++found_neghbours_count_;
      }
      ++current_cell;
  }
  clock_t after_traversal = clock();
  perf_tester_.full_traversal_time += after_traversal - before_traversal;
  clock_t finish = clock();
  perf_tester_.full_search_time += finish - start;
}

template<>
inline void RecordToMetainfoAndDistance<RerankADC8, PointId>(const Coord* point, const RerankADC8& record,
                                                             pair<Distance, PointId>* result,
                                                             const vector<int>& cell_coordinates,
                                                             const vector<Centroids>& fine_vocabs) {
  result->second = record.pid;
  int coarse_clusters_count = cell_coordinates.size();
  int fine_clusters_count = fine_vocabs.size();
  int coarse_to_fine_ratio = fine_clusters_count / coarse_clusters_count;
  int subvectors_dim = SPACE_DIMENSION / fine_clusters_count;
  char* rerank_info_ptr = (char*)&record + sizeof(record.pid);
  for(int centroid_index = 0; centroid_index < fine_clusters_count; ++centroid_index) {
    int start_dim = centroid_index * subvectors_dim;
    int final_dim = start_dim + subvectors_dim;
    FineClusterId pid_nearest_centroid = *((FineClusterId*)rerank_info_ptr);
    rerank_info_ptr += sizeof(FineClusterId);
    int current_coarse_index = centroid_index / coarse_to_fine_ratio;
    Distance subvector_distance = 0;
    for(int i = start_dim; i < final_dim; ++i) {
      Coord diff = fine_vocabs[centroid_index][pid_nearest_centroid][i - start_dim] - point[i];
        subvector_distance += diff * diff;
    }
    result->first += subvector_distance;
  }
}

template<>
inline void RecordToMetainfoAndDistance<RerankADC16, PointId>(const Coord* point, const RerankADC16& record,
                                                              pair<Distance, PointId>* result,
                                                              const vector<int>& cell_coordinates,
                                                              const vector<Centroids>& fine_vocabs) {
  result->second = record.pid;
  int coarse_clusters_count = cell_coordinates.size();
  int fine_clusters_count = fine_vocabs.size();
  int coarse_to_fine_ratio = fine_clusters_count / coarse_clusters_count;
  int subvectors_dim = SPACE_DIMENSION / fine_clusters_count;
  char* rerank_info_ptr = (char*)&record + sizeof(record.pid);
  for(int centroid_index = 0; centroid_index < fine_clusters_count; ++centroid_index) {
    int start_dim = centroid_index * subvectors_dim;
    int final_dim = start_dim + subvectors_dim;
    FineClusterId pid_nearest_centroid = *((FineClusterId*)rerank_info_ptr);
    rerank_info_ptr += sizeof(FineClusterId);
    int current_coarse_index = centroid_index / coarse_to_fine_ratio;
    Distance subvector_distance = 0;
    for(int i = start_dim; i < final_dim; ++i) {
      Coord diff = fine_vocabs[centroid_index][pid_nearest_centroid][i - start_dim] - point[i];
      subvector_distance += diff * diff;
    }
    result->first += subvector_distance;
  }
}

template class Searcher<RerankADC8, PointId>;
template class Searcher<RerankADC16, PointId>;
//template class Searcher<PointId, PointId>;

#endif

