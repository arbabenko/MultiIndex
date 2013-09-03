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
  //PerfTester& GetPerfTester();
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
  //void InitBlasStructures();
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
  ifstream multi_array(string(index_files_prefix + "_multi_array.bin").c_str(), ios::binary);
  if(!multi_array.good()) {
    throw std::logic_error("Bad input multiarray stream");
  }
  boost::archive::binary_iarchive arc_multi_array(multi_array);
  arc_multi_array >> multiindex_.multiindex;
  for (int i=0;i<16384;++i) cout << multiindex_.cell_edges.table[i] << endl;
  cout << "Multiindex deserialized...\n";
  ifstream norms(string(index_files_prefix + "_prec_norms.bin").c_str(), ios::binary);
  if(!norms.good()) {
    throw std::logic_error("Bad input norms stream");
  }
  boost::archive::binary_iarchive arc_norms(norms);
  arc_norms >> precomputed_norms_;
  cout << "Multiindex deserialized...\n";

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
  //merger_.GetYieldedItems().table.resize(std::pow((float)subspace_centroids_to_consider,
	//	                                         (int)coarse_vocabs_.size()));
  //for(int i = 0; i < coarse_vocabs_.size(); ++i) {
  //  merger_.GetYieldedItems().dimensions.push_back(subspace_centroids_to_consider);
  //}
  //InitBlasStructures();
}

//template<class Record, class MetaInfo>
//void MultiSearcher<Record, MetaInfo>::InitBlasStructures(){
//  coarse_vocabs_matrices_.resize(coarse_vocabs_.size());
//  coarse_centroids_norms_.resize(coarse_vocabs_.size(), vector<float>(coarse_vocabs_[0].size()));
//  for(int coarse_id = 0; coarse_id < coarse_vocabs_matrices_.size(); ++coarse_id) {
//    coarse_vocabs_matrices_[coarse_id] = new float[coarse_vocabs_[0].size() * coarse_vocabs_[0][0].size()];
//    for(int i = 0; i < coarse_vocabs_[0].size(); ++i) {
//      Coord norm = 0;
//      for(int j = 0; j < coarse_vocabs_[0][0].size(); ++j) {
//        coarse_vocabs_matrices_[coarse_id][coarse_vocabs_[0][0].size() * i + j] = coarse_vocabs_[coarse_id][i][j];
//        norm += coarse_vocabs_[coarse_id][i][j] * coarse_vocabs_[coarse_id][i][j];
//      }
//      coarse_centroids_norms_[coarse_id][i] = norm;
//    }
//  }
//  products_ = new Coord[coarse_vocabs_[0].size()];
//  query_norms_.resize(coarse_vocabs_[0].size());
//  residual_ = new Coord[coarse_vocabs_[0][0].size() * coarse_vocabs_.size()];
//}

//template<class Record, class MetaInfo>
//PerfTester& MultiSearcher<Record, MetaInfo>::GetPerfTester() {
//  return perf_tester_;
//}

//template<class Record, class MetaInfo>
//void MultiSearcher<Record, MetaInfo>::GetNearestSubspacesCentroids(const Point& point,
//                                                                   const int subspace_centroins_count,
//                                                                   vector<NearestSubspaceCentroids>*
//                                                                   subspaces_short_lists) const {
//  std::stringstream aa;
//  subspaces_short_lists->resize(coarse_vocabs_.size());
//  Dimensions subspace_dimension = point.size() / coarse_vocabs_.size();
//  for(int subspace_index = 0; subspace_index < coarse_vocabs_.size(); ++subspace_index) {
//    Dimensions start_dim = subspace_index * subspace_dimension;
//    Dimensions final_dim = std::min((Dimensions)point.size(), start_dim + subspace_dimension);
//    Coord query_norm = cblas_sdot(final_dim - start_dim, &(point[start_dim]), 1, &(point[start_dim]), 1);
//    std::fill(query_norms_.begin(), query_norms_.end(), query_norm);
//    cblas_saxpy(coarse_vocabs_[0].size(), 1, &(coarse_centroids_norms_[subspace_index][0]), 1, &(query_norms_[0]), 1);
//    cblas_sgemv(CblasRowMajor, CblasNoTrans, coarse_vocabs_[0].size(), subspace_dimension, -2.0,
//                coarse_vocabs_matrices_[subspace_index], subspace_dimension, &(point[start_dim]), 1, 1, &(query_norms_[0]), 1);
//    subspaces_short_lists->at(subspace_index).resize(query_norms_.size());
//    for(int i = 0; i < query_norms_.size(); ++i) {
//      subspaces_short_lists->at(subspace_index)[i] = std::make_pair(query_norms_[i], i);
//    }
//    std::nth_element(subspaces_short_lists->at(subspace_index).begin(),
//                     subspaces_short_lists->at(subspace_index).begin() + subspace_centroins_count,
//                     subspaces_short_lists->at(subspace_index).end());
//    subspaces_short_lists->at(subspace_index).resize(subspace_centroins_count);
//    std::sort(subspaces_short_lists->at(subspace_index).begin(),
//              subspaces_short_lists->at(subspace_index).end());
//  }
//}

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

//template<class Record, class MetaInfo>
//bool MultiSearcher<Record, MetaInfo>::TraverseNextMultiIndexCell(const Point& point,
//                                                                 vector<pair<Distance, MetaInfo> >*
//                                                                             nearest_subpoints) const {
//  MergedItemIndices cell_inner_indices;
//  clock_t before = clock();
//  if(!merger_.GetNextMergedItemIndices(&cell_inner_indices)) {
//    return false;
//  }
//  clock_t after = clock();
//  perf_tester_.cell_coordinates_time += after - before;
//  vector<int> cell_coordinates(cell_inner_indices.size());
//  for(int list_index = 0; list_index < merger_.lists_ptr->size(); ++list_index) {
//    cell_coordinates[list_index] = merger_.lists_ptr->at(list_index)[cell_inner_indices[list_index]].second;
//  }
//  int cell_start, cell_finish;
//  before = clock();
//  GetCellEdgesInMultiIndexArray(cell_coordinates, &cell_start, &cell_finish);
//  after = clock();
//  perf_tester_.cell_edges_time += after - before;
//  if(cell_start >= cell_finish) {
//    return true;
//  }
//  typename vector<Record>::const_iterator it = multiindex_.multiindex.begin() + cell_start;
//  GetResidual(point, cell_coordinates, coarse_vocabs_, residual_);
//  cell_finish = std::min((int)cell_finish, cell_start + (int)nearest_subpoints->size() - found_neghbours_count_);
//  for(int array_index = cell_start; array_index < cell_finish; ++array_index) {
//    if(rerank_mode_ == USE_RESIDUALS) {
//      RecordToMetainfoAndDistance<Record, MetaInfo>(residual_, *it,
//                                                    &(nearest_subpoints->at(found_neghbours_count_)),
//                                                    cell_coordinates, fine_vocabs_);
//    } else if(rerank_mode_ == USE_INIT_POINTS) {
//      RecordToMetainfoAndDistance<Record, MetaInfo>(&(point[0]), *it,
//                                                    &(nearest_subpoints->at(found_neghbours_count_)),
//                                                    cell_coordinates, fine_vocabs_);
//    }
//    perf_tester_.NextNeighbour();
//    ++found_neghbours_count_;
//    ++it;
//  }
//  return true;
//}


template<class Record, class MetaInfo>
void Searcher<Record, MetaInfo>::GetNearestNeighbours(const Point& point, int k, 
                                                           vector<pair<Distance, MetaInfo> >* neighbours) const {
  found_neghbours_count_ = 0;
  neighbours->resize(k);
  vector<float> main_products(main_vocabs_.size());
  vector<float> res_products(res_vocabs_.size());
  vector<pair<float, int> > temp(main_products.size());
  for (int i = 0; i < main_vocabs_.size(); ++i) {
    main_products[i] = cblas_sdot(point.size(), &(point[0]), 1, &(main_vocabs_[i][0]), 1);
    float norm = cblas_sdot(point.size(), &(main_vocabs_[i][0]), 1, &(main_vocabs_[i][0]), 1);
    temp[i] = std::make_pair(norm / 2 - main_products[i], i);
  }
  
  for (int i = 0; i < res_vocabs_.size(); ++i) {
    res_products[i] = cblas_sdot(point.size(), &(point[0]), 1, &(res_vocabs_[i][0]), 1);
  }
  std::sort(temp.begin(), temp.end());
  std::multimap<float, pair<ClusterId, ClusterId> > queue;
  for(int i = 0; i < main_centroids_to_consider_; ++i) {
    for(int j = 0; j < res_vocabs_.size(); ++j) {
      int main_cluster = temp[i].second;
      float score = precomputed_norms_[main_cluster][j] - 2 * main_products[main_cluster] - 2 * res_products[j];
      queue.insert(std::make_pair(score, std::make_pair(main_cluster,j)));
    }
  }
  cout << "Traverse started" << endl;
  std::multimap<float, pair<ClusterId, ClusterId> >::iterator current_cell = queue.begin();
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
  //assert(k > 0);
  //perf_tester_.handled_queries_count += 1;
  //neighbours->resize(k);
  //perf_tester_.ResetQuerywiseStatistic();
  //clock_t start = clock();
  //perf_tester_.search_start = start;
  //clock_t before = clock();
  //vector<NearestSubspaceCentroids> subspaces_short_lists;
  //assert(subspace_centroids_to_consider_ > 0);
  //GetNearestSubspacesCentroids(point, subspace_centroids_to_consider_, &subspaces_short_lists);
  //clock_t after = clock();
  //perf_tester_.nearest_subcentroids_time += after - before;
  //clock_t before_merger = clock();
  //merger_.setLists(subspaces_short_lists);
  //clock_t after_merger = clock();
  //perf_tester_.merger_init_time += after_merger - before_merger;
  //clock_t before_traversal = clock();
  //found_neghbours_count_ = 0;
  //bool traverse_next_cell = true;
  //int cells_visited = 0;
  //while(found_neghbours_count_ < k && traverse_next_cell) {
  //  perf_tester_.cells_traversed += 1;
  //  traverse_next_cell = TraverseNextMultiIndexCell(point, neighbours);
  //  cells_visited += 1;
  //}
  //clock_t after_traversal = clock();
  //perf_tester_.full_traversal_time += after_traversal - before_traversal;
  //if(do_rerank_) {
  //  //if(neighbours->size() > 10000) {
  //  //  std::nth_element(neighbours->begin(), neighbours->begin() + 10000, neighbours->end());
  //  //  neighbours->resize(10000);
  //  //}
  //  //std::sort(neighbours->begin(), neighbours->end());
  //}
  //clock_t finish = clock();
  //perf_tester_.full_search_time += finish - start;
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

