/** @file */

// Copyright 2012 Yandex Artem Babenko
#ifndef INDEXER_H_
#define INDEXER_H_

#include <ctime>
#include <limits>
#include <map>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>

#include <boost/lexical_cast.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>

#include "data_util.h"
#include "multitable.hpp"


using std::ifstream;
using std::map;
using std::multimap;
using std::ofstream;
using std::string;

using boost::lexical_cast;
using boost::split;

extern int THREADS_COUNT;

extern Dimensions SPACE_DIMENSION;

extern enum PointType point_type;

IndexConfig gConfig;

/**
 * This is the main class for creating hierarchical index for a set of points
 * in a multidimensional space. Clusterization and vocabs learning are maid
 * outside of this class, indexer receives prepared vocabs in input
 */
template<class Record>
class Indexer {
 public:
 /**
  * This is the simple Indexer constructor
  */
  Indexer();
 /**
  * This is the main function of Indexer
  * @param points_filename file with points in .fvecs or .bvecs format
  * @param points_count how many points should we index
  * @param main_vocabs vocabularies for initial points
  * @param res_vocabs vocabularies for residuals
  * @param mode determines the way of rerank info calculating
  * @param build_quantization should we get quantization or not
  * @param files_prefix all index filenames will have this prefix
  * @param quantization_filename file with quantization (if exists)
  */
  void BuildHierIndex(const string& points_filename,
                      const string& metainfo_filename,
                      const int points_count,
                      const Centroids& main_vocabs,
                      const Centroids& res_vocabs,
                      const RerankMode& mode,
                      const bool build_quantization,
                      const string& files_prefix,
                      const string& quantization_filename,
                      int main_centroids_count);
 private:
 /**
  * This function calculates for each point its main vocab and residual vocab
  * @param points_filename file with points in .fvecs or .bvecs format
  * @param points_count how many points should we handle
  * @param main_vocabs vocabularies for initial points
  * @param res_vocabs vocabularies for residuals
  */
  void PrepareQuantization(const string& points_filename,
                           const int points_count,
                           const Centroids& main_vocabs,
                           const Centroids& res_vocabs);
 /**
  * This function prepares for each point in subset its main vocab and residual vocab
  * @param points_filename file with points in .fvecs or .bvecs format
  * @param start_pid identifier of the first point in subset
  * @param subset_size points count in subset
  * @param main_vocabs vocabularies for initial points
  * @param res_vocabs vocabularies for residuals
  * @param transposed_quantizations result
  */
  void GetQuantizationsForSubset(const string& points_filename,
                                 const int start_pid,
                                 const int subset_size,
                                 const Centroids& main_vocabs,
                                 const Centroids& res_vocabs,
                                 vector<vector<ClusterId> >*
                                 transposed_quantizations);
 /**
  * This function serializes prepared quantizations to file
  * @param transposed_quantizations quantizations to serialize.
  * They are transposed because of effective memory usage
  * @param filename file we should serialize to
  */
  void SerializeQuantizations(const vector<vector<ClusterId> >&
                              transposed_quantizations,
                              const string& filename);
 /**
  * This function saves hierarchical index to files.
  * All filenames start form the common files prefix
  */
  void SerializeHierIndexFiles();
 /**
  * This function converts counts of points in cells to cell edges
  */
  void ConvertPointsInCellsCountToCellEdges();

 /**
  * This function fills hierarchical index data structures.
  * @param points_filename file with points in .fvecs or .bvecs format
  * @param points_count how many points should we index
  * @param main_vocabs vocabularies for initial points
  * @param res_vocabs vocabularies for residuals
  * @param mode determines the way of rerank info calculating
  */
  void FillHierIndex(const string& points_filename,
                     const int points_count,
                     const Centroids& main_vocabs,
                     const Centroids& res_vocabs,
                     const RerankMode& mode);
 /**
  * This function fills hierarchical index data structures.
  * @param points_filename file with points in .fvecs or .bvecs format
  * @param start_pid identifier of the first point in subset
  * @param subset_size points count in subset
  * @param main_vocabs vocabularies for initial points
  * @param res_vocabs vocabularies for residuals
  * @param mode determines the way of rerank info calculating
  * @param points_written_in_index auxillary structure for correct index filling
  */
  void FillHierIndexForSubset(const string& points_filename,
                              const PointId start_pid,
                              const int points_count,
                              const Centroids& main_vocabs,
                              const Centroids& res_vocabs,
                              const RerankMode& mode,
                              Multitable<int>* points_written_in_index);

 /**
  * This function reads point quantization from file
  * @param pid identifier of target point
  * @param filename file with quantizations
  * @param quantization result
  */
  void GetPointQuantization(const PointId pid, const string& filename,
                            vector<ClusterId>* quantization);
 /**
  * This function restores counts of points from their quantizations
  * @param points_filename file with points in .fvecs or .bvecs format
  * @param points_count how many points should we index
  * @param main_vocabs vocabularies for initial points
  * @param res_vocabs vocabularies for residuals
  * We need them to init counts table correctly
  */
  void RestorePointsInCellsCountFromQuantization(const string& points_filename,
                                                 const int points_count,
                                                 const Centroids& main_vocabs,
                                                 const Centroids& res_vocabs);

  void PrecomputeEffectiveCentroidsNorms(const Centroids& main_vocabs,
                                         const Centroids& res_vocabs);
 /**
  * This simple function returns size of one coordinate of input point
  */
  int GetInputCoordSizeof();
 /**
  * This simple function reads one point from input stream
  * @param input input stream
  * @param point result point
  */
  void ReadPoint(ifstream& input, Point* point);

  void InitBlasStructures(const Centroids& main_vocabs, const Centroids& res_vocabs);

 /**
  *  All index filenames will start from this prefix
  */
  string files_prefix_;
 /**
  *  Filename of file with quantizations
  */
  string quantization_filename_;
 /**
  *  Table with number of points in each cell
  */
  Multitable<int> point_in_cells_count_;
 /**
  *  Multiindex
  */
  MultiIndex<Record> multiindex_;
 /**
  *  Mutex for critical section in filling index stage
  */
  boost::mutex cell_counts_mutex_;

  vector<float> main_norms_;

  vector<vector<float> > precomputed_norms_;

  int main_centroids_count_;

  float* main_vocabs_matrix_;

  float* res_vocabs_matrix_;
};

template<class Record>
inline void GetRecord(const Point& point, const PointId pid,
                      const vector<ClusterId> quantization,
                      const Centroids& main_vocabs,
                      const Centroids& res_vocabs,
                      Record* result) {
}

//////////////////// IMPLEMENTATION //////////////////////
template<class Record>
Indexer<Record>::Indexer() {
}

template<class Record>
int Indexer<Record>::GetInputCoordSizeof() {
  if(point_type == FVEC) {
    return (int)sizeof(float);
  } else if(point_type == BVEC) {
    return (int)sizeof(unsigned char);
  }
}

template<class Record>
void Indexer<Record>::ReadPoint(ifstream& input, Point* point) {
  if(!input.good()) {
    throw std::logic_error("Bad input stream");
  }
  if(point_type == FVEC) {
    ReadVector<float, Coord>(input, point);
  } else if(point_type == BVEC) {
    ReadVector<unsigned char, Coord>(input, point);
  }    
}

template<class Record>
void Indexer<Record>::PrecomputeEffectiveCentroidsNorms(const Centroids& main_vocabs,
                                                        const Centroids& res_vocabs) {
  main_norms_.resize(main_vocabs.size());
  for(int i = 0; i < main_vocabs.size(); ++i) {
    main_norms_[i] = cblas_sdot(main_vocabs[i].size(), &(main_vocabs[i][0]), 1, &(main_vocabs[i][0]), 1);
  }
  precomputed_norms_.resize(main_vocabs.size());
  for(int i = 0; i < main_vocabs.size(); ++i) {
    precomputed_norms_[i].resize(res_vocabs.size());
    for(int j = 0; j < res_vocabs.size(); ++j) {
      Point temp = main_vocabs[i];
      cblas_saxpy(temp.size(), 1, &(res_vocabs[j][0]), 1, &(temp[0]), 1);
      precomputed_norms_[i][j] = cblas_sdot(temp.size(), &(temp[0]), 1, &(temp[0]), 1);
    }
  }
}


template<class Record>
void Indexer<Record>::SerializeQuantizations(const vector<vector<ClusterId> >&
		                                         transposed_quantizations,
                                             const string& filename) {
  ofstream quantizations_stream;
  quantizations_stream.open(filename.c_str(), ios::binary);
  if(!quantizations_stream.good()) {
    throw std::logic_error("Bad input stream");
  }
  cout << "Writing quantizations started " << filename << endl;
  for(PointId pid = 0; pid < transposed_quantizations[0].size(); ++pid) {
    for(int qauntization_idx = 0;
        qauntization_idx < transposed_quantizations.size();
        ++qauntization_idx) {
      ClusterId quantization = transposed_quantizations[qauntization_idx][pid];
      quantizations_stream.write((char*)&quantization, sizeof(quantization));
    }
  }
  quantizations_stream.close();
  cout << "Writing quantizations started" << endl;
}

template<class Record>
void Indexer<Record>::SerializeHierIndexFiles() {
  cout << "Start hierindex serializing....\n";
  ofstream cell_edges(string(files_prefix_ + "_cell_edges.bin").c_str(), ios::binary);
  boost::archive::binary_oarchive arc_cell_edges(cell_edges);
  arc_cell_edges << multiindex_.cell_edges;
  ofstream multi_array(string(files_prefix_ + "_multi_array.bin").c_str(), ios::binary);
  boost::archive::binary_oarchive arc_multi_array(multi_array);
  arc_multi_array << multiindex_.multiindex;
  ofstream norms(string(files_prefix_ + "_prec_norms.bin").c_str(), ios::binary);
  boost::archive::binary_oarchive arc_norms(norms);
  arc_norms << precomputed_norms_;
  ofstream main_norms(string(files_prefix_ + "_main_norms.bin").c_str(), ios::binary);
  boost::archive::binary_oarchive arc_main_norms(main_norms);
  arc_main_norms << main_norms_;
  cout << "Finish hierindex serializing....\n";
}

//std::ofstream dist("distor_hier.txt");

template<class Record>
void Indexer<Record>::GetQuantizationsForSubset(const string& points_filename,
                                                const int start_pid,
                                                const int subset_size,
                                                const Centroids& main_vocabs,
                                                const Centroids& res_vocabs,
                                                vector<vector<ClusterId> >*
                                                      transposed_quantizations) {
  ifstream point_stream;
  point_stream.open(points_filename.c_str(), ios::binary);
  if(!point_stream.good()) {
    throw std::logic_error("Bad input points stream");
  }
  // we assume points are stored in .fvecs or .bvecs format
  point_stream.seekg(start_pid * (GetInputCoordSizeof() * SPACE_DIMENSION + sizeof(Dimensions)), ios::beg);
  vector<ClusterId> quantization(2);
  for(int point_number = 0; point_number < subset_size; ++point_number) {
    if(point_number % 10000 == 0) {
      cout << "Getting quantization, point # " << start_pid + point_number << endl;
    }
    Point current_point;
    ReadPoint(point_stream, &current_point);
    ////// GETTING QUANTIZATIONS
    vector<float> main_products(main_vocabs.size(), 0);
    vector<float> res_products(res_vocabs.size(), 0);
    vector<pair<float, int> > distance_to_clusterId(main_products.size());
    cblas_sgemv(CblasRowMajor, CblasNoTrans, main_vocabs.size(), main_vocabs[0].size(), 1.0,
                main_vocabs_matrix_, main_vocabs[0].size(), &(current_point[0]), 1, 1, &(main_products[0]), 1);
    for (int main_cid = 0; main_cid < main_vocabs.size(); ++main_cid) {
      //main_products[main_cid] = cblas_sdot(current_point.size(), &(current_point[0]),
      //                                     1, &(main_vocabs[main_cid][0]), 1);
      distance_to_clusterId[main_cid] = std::make_pair(main_norms_[main_cid] / 2 - main_products[main_cid], main_cid);
    }
    //for (int res_cid = 0; res_cid < res_vocabs.size(); ++res_cid) {
    //  res_products[res_cid] = cblas_sdot(current_point.size(), &(current_point[0]),
    //                                       1, &(res_vocabs[res_cid][0]), 1);
    //}
    cblas_sgemv(CblasRowMajor, CblasNoTrans, res_vocabs.size(), res_vocabs[0].size(), 1.0,
                res_vocabs_matrix_, res_vocabs[0].size(), &(current_point[0]), 1, 1, &(res_products[0]), 1);
    std::sort(distance_to_clusterId.begin(), distance_to_clusterId.end());

    ClusterId optimal_main = 0;
    ClusterId optimal_res = 0;

    float optimal_distance = std::numeric_limits<float>::infinity();
    for(int main_cid = 0; main_cid < main_centroids_count_; ++main_cid) {
      for(int res_cid = 0; res_cid < res_vocabs.size(); ++res_cid) {
        int current_optimal_main_cid = distance_to_clusterId[main_cid].second;
        float distance = precomputed_norms_[current_optimal_main_cid][res_cid] -
                         2 * main_products[current_optimal_main_cid] -
                         2 * res_products[res_cid];
        if(distance < optimal_distance) {
          optimal_main = current_optimal_main_cid;
          optimal_res = res_cid;
          optimal_distance = distance;
        }
      }
    }
    transposed_quantizations->at(0)[start_pid + point_number] = optimal_main;
    quantization[0] = optimal_main;
    transposed_quantizations->at(1)[start_pid + point_number] = optimal_res;
    quantization[1] = optimal_res;
    // just for testing
    //cblas_saxpy(current_point.size(), -1, &(main_vocabs[optimal_main][0]), 1, &(current_point[0]), 1);
    //float distors1 = cblas_sdot(current_point.size(), &(current_point[0]), 1, &(current_point[0]), 1);
    //cblas_saxpy(current_point.size(), -1, &(res_vocabs[optimal_res][0]), 1, &(current_point[0]), 1);
    //float distors2 = cblas_sdot(current_point.size(), &(current_point[0]), 1, &(current_point[0]), 1);
    //
    int global_index = point_in_cells_count_.GetCellGlobalIndex(quantization);
    cell_counts_mutex_.lock();
    //dist << distors1 << " " << distors2 << std::endl; // debugging
    ++(point_in_cells_count_.table[global_index]);
    cell_counts_mutex_.unlock();
  }
}

template<class Record>
void Indexer<Record>::PrepareQuantization(const string& points_filename,
                                          const int points_count,
                                          const Centroids& main_vocabs,
                                          const Centroids& res_vocabs) {
  // we use transposed quantizations for efficient memory usage
  vector<vector<ClusterId> > transposed_quantizations; 
  transposed_quantizations.resize(2);
  transposed_quantizations[0].resize(points_count);
  transposed_quantizations[1].resize(points_count);
  vector<int> multiindex_table_dimensions;
  multiindex_table_dimensions.push_back(main_vocabs.size());
  multiindex_table_dimensions.push_back(res_vocabs.size());
  point_in_cells_count_.Resize(multiindex_table_dimensions);
  cout << "Memory for quantizations allocated" << endl;
  boost::thread_group index_threads;
  int thread_points_count = points_count / THREADS_COUNT;
  for(int thread_id = 0; thread_id < THREADS_COUNT; ++thread_id) {
    PointId start_pid = thread_points_count * thread_id;
    index_threads.create_thread(boost::bind(&Indexer::GetQuantizationsForSubset,
                                            this, points_filename, start_pid, thread_points_count,
                                            boost::cref(main_vocabs), boost::cref(res_vocabs),
                                            &transposed_quantizations));
  }
  index_threads.join_all();
  if(quantization_filename_.empty()) {
    quantization_filename_ = files_prefix_ + "_quantizations.bin";
  }
  cout << "Quantizations are calculated" << endl;
  SerializeQuantizations(transposed_quantizations, quantization_filename_);
  cout << "Quantizations are serialized" << endl;
}

template<class Record>
void Indexer<Record>::ConvertPointsInCellsCountToCellEdges() {
  cout << "Converting points in cells to cell edges...\n";
  multiindex_.cell_edges = point_in_cells_count_;
  multiindex_.cell_edges.table[0] = 0;
  for(int global_index = 1;
      global_index < point_in_cells_count_.table.size();
      ++global_index) {
    multiindex_.cell_edges.table[global_index] = multiindex_.cell_edges.table[global_index - 1] +
                                                 point_in_cells_count_.table[global_index - 1];
  }
  // we do not need this table more
  point_in_cells_count_.table.clear();
  cout << "Finish converting points in cells to cell edges...\n";
}

template<class Record>
void Indexer<Record>::GetPointQuantization(const PointId pid,
                                           const string& filename,
                                           vector<ClusterId>* quantization) {
  ifstream quantization_stream;
  quantization_stream.open(filename.c_str(), ios::binary);
  if(!quantization_stream.good()) {
    throw std::logic_error("Bad input quantizations stream");
  }
  quantization_stream.seekg((long long)pid * sizeof(ClusterId) * 2, ios::beg);
  quantization_stream.read((char*)&(quantization->at(0)), sizeof(quantization->at(0)));
  quantization_stream.read((char*)&(quantization->at(1)), sizeof(quantization->at(1)));
}

template<class Record>
void Indexer<Record>::FillHierIndexForSubset(const string& points_filename,
                                             const PointId start_pid,
                                             const int points_count,
                                             const Centroids& main_vocabs,
                                             const Centroids& res_vocabs,
                                             const RerankMode& mode,
                                             Multitable<int>* points_written_in_index) {
  ifstream point_stream;
  point_stream.open(points_filename.c_str(), ios::binary);
  if(!point_stream.good()) {
    throw std::logic_error("Bad input points stream");
  }
  point_stream.seekg((long long)start_pid * (GetInputCoordSizeof() * SPACE_DIMENSION + sizeof(Dimensions)), ios::beg);
  for(int point_number = 0; point_number < points_count; ++point_number) {
    if(point_number % 10000 == 0) {
      cout << "Filling hierindex, point # " << start_pid + point_number << endl;
    }
    Point current_point;
    ReadPoint(point_stream, &current_point);
    vector<ClusterId> quantization(2);
    GetPointQuantization(start_pid + point_number, quantization_filename_, &quantization);

    cell_counts_mutex_.lock();
    int current_written_count = points_written_in_index->GetValue(quantization);
    int pid_multiindex = multiindex_.cell_edges.GetValue(quantization) + current_written_count;
    multiindex_.multiindex[pid_multiindex].pid = start_pid + point_number;
    // TODO - add rerank info
    //GetRecord<Record>(current_point, start_pid + point_number,
    //                  coarse_quantization, coarse_vocabs, &(multiindex_.multiindex[pid_multiindex]));
    points_written_in_index->SetValue(current_written_count + 1, quantization);
    cell_counts_mutex_.unlock();
  }
}

template<class Record>
void Indexer<Record>::FillHierIndex(const string& points_filename,
                                    const int points_count,
                                    const Centroids& main_vocabs,
                                    const Centroids& res_vocabs,
                                    const RerankMode& mode) {
  ConvertPointsInCellsCountToCellEdges();
  multiindex_.multiindex.resize(points_count);
  cout << "Indexing started..." << endl;

  Multitable<int> points_written_in_index(multiindex_.cell_edges.dimensions);
  int thread_points_count = points_count / THREADS_COUNT;
  boost::thread_group threads;
  for(int thread_id = 0; thread_id < THREADS_COUNT; ++thread_id) {
    PointId start_pid = thread_points_count * thread_id;
    threads.create_thread(boost::bind(&Indexer::FillHierIndexForSubset, this, points_filename, start_pid,
                                      thread_points_count, boost::cref(main_vocabs),
                                      boost::cref(res_vocabs), mode, &points_written_in_index));
  }
  threads.join_all();
  cout << "Indexing finished..." << endl;
}

template<class Record>
void Indexer<Record>::RestorePointsInCellsCountFromQuantization(const string& points_filename,
                                                                const int points_count,
                                                                const Centroids& main_vocabs,
                                                                const Centroids& res_vocabs) {
  vector<int> dimensions;
  dimensions.push_back(main_vocabs.size());
  dimensions.push_back(res_vocabs.size());
  point_in_cells_count_.Resize(dimensions);
  ifstream quantization_stream;
  quantization_stream.open(quantization_filename_.c_str(), ios::binary);
  if(!quantization_stream.good()) {
    throw std::logic_error("Bad input quantizations stream");
  }
  vector<ClusterId> quantization(2);
  for(PointId pid = 0; pid < points_count; ++pid) {
    if(pid % 100000 == 0) {
      cout << pid << endl;
    }
    quantization_stream.read((char*)&(quantization[0]), sizeof(ClusterId));
    quantization_stream.read((char*)&(quantization[1]), sizeof(ClusterId));
    cout << quantization[0] << " " << quantization[1] << endl;
    int cell_global_index = point_in_cells_count_.GetCellGlobalIndex(quantization);
    point_in_cells_count_.table[cell_global_index] += 1;
  }
}

template<class Record>
void Indexer<Record>::BuildHierIndex(const string& points_filename,
                                     const string& metainfo_filename,
                                     const int points_count,
                                     const Centroids& main_vocabs,
                                     const Centroids& res_vocabs,
                                     const RerankMode& mode,
                                     const bool build_quantization,
                                     const string& files_prefix,
                                     const string& quantization_filename,
                                     int main_centroids_count) {
  //InitParameters<Record>(fine_vocabs, mode, metainfo_filename);
  InitBlasStructures(main_vocabs, res_vocabs);
  files_prefix_ = files_prefix;
  main_centroids_count_ = main_centroids_count;
  quantization_filename_ = quantization_filename;
  PrecomputeEffectiveCentroidsNorms(main_vocabs, res_vocabs);
  cout << "Norms are precomputed" << endl;
  if(build_quantization) {
    PrepareQuantization(points_filename, points_count, main_vocabs, res_vocabs);
  } else {
  RestorePointsInCellsCountFromQuantization(points_filename,
                                            points_count,
                                            main_vocabs,
                                            res_vocabs);
  }
  FillHierIndex(points_filename, points_count, main_vocabs, res_vocabs, mode);
  cout << "Hierindex created" << endl;
  SerializeHierIndexFiles();
  cout << "Hierindex serialized" << endl;
}

template<class Record>
void Indexer<Record>::InitBlasStructures(const Centroids& main_vocabs,
                                         const Centroids& res_vocabs) {
  main_vocabs_matrix_ = new float[main_vocabs.size() * main_vocabs[0].size()];
  for(int cid = 0; cid < main_vocabs.size(); ++cid) {
    for(int coord = 0; coord < main_vocabs[0].size(); ++coord) {
        main_vocabs_matrix_[cid * main_vocabs[0].size() + coord] = main_vocabs[cid][coord];
    }
  }
  res_vocabs_matrix_ = new float[res_vocabs.size() * res_vocabs[0].size()];
  for(int cid = 0; cid < res_vocabs.size(); ++cid) {
    for(int coord = 0; coord < res_vocabs[0].size(); ++coord) {
        res_vocabs_matrix_[cid * res_vocabs[0].size() + coord] = res_vocabs[cid][coord];
    }
  }
}

template<>
inline void GetRecord<PointId> (const Point& point, const PointId pid,
                                const vector<ClusterId> quantization,
                                const Centroids& coarse_vocabs,
                                const Centroids& res_vocabs,
                                PointId* result) {
  *result = pid;
}

inline void FillAdcInfo(const Point& point, const PointId pid,
                        const vector<Centroids>& rerank_vocabs,
                        char* result) {
  int subvectors_count = rerank_vocabs.size();
  int subvector_dim = point.size() / subvectors_count;
  for(int subvector_index = 0; subvector_index < subvectors_count; ++subvector_index) {
    Dimensions start_dim = subvector_index * subvector_dim;
    Dimensions final_dim = start_dim + subvector_dim;
    *((FineClusterId*)result) = (FineClusterId)GetNearestClusterId(point, rerank_vocabs[subvector_index],
			                                                             start_dim, final_dim);
    result += sizeof(FineClusterId);
  }
}

template<>
inline void GetRecord<RerankADC8>(const Point& point, const PointId pid,
                                  const vector<ClusterId> quantization,
                                  const Centroids& coarse_vocabs,
                                  const Centroids& res_vocabs,
                                  RerankADC8* result) {
  result->pid = pid;
  char* rerank_info_ptr = (char*)result + sizeof(pid);
  if(gConfig.rerank_mode == USE_RESIDUALS) {
    Point residual;
    //GetResidual(point, quantization, coarse_vocabs, &residual);
    FillAdcInfo(residual, pid, gConfig.fine_vocabs, rerank_info_ptr);
  } else if (gConfig.rerank_mode == USE_INIT_POINTS) {
    FillAdcInfo(point, pid, gConfig.fine_vocabs, rerank_info_ptr);
  }
}

template<>
inline void GetRecord<RerankADC16> (const Point& point, const PointId pid,
                                    const vector<ClusterId> quantization,
                                    const Centroids& coarse_vocabs,
                                    const Centroids& res_vocabs,
                                    RerankADC16* result) {
  result->pid = pid;
  char* rerank_info_ptr = (char*)result + sizeof(pid);
  if(gConfig.rerank_mode == USE_RESIDUALS) {
    Point residual;
    //GetResidual(point, quantization, coarse_vocabs, &residual);
    FillAdcInfo(residual, pid, gConfig.fine_vocabs, rerank_info_ptr);
  } else if (gConfig.rerank_mode == USE_INIT_POINTS) {
    FillAdcInfo(point, pid, gConfig.fine_vocabs, rerank_info_ptr);
  }
}

#endif


