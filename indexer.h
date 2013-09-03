/** @file */

// Copyright 2012 Yandex Artem Babenko
#ifndef INDEXER_H_
#define INDEXER_H_

#include <ctime>
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
 * in a multidimensional space. Clusterization and vocabs learning happen
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
  * @param build_coarse_quantization should we get coarse quantization or not
  * @param files_prefix all index filenames will have this prefix
  * @param coarse_quantization_filename file with coarse quantization (if exists)
  */
  void BuildHierIndex(const string& points_filename,
                      const string& metainfo_filename,
                      const int points_count,
                      const Centroids& main_vocabs,
                      const Centroids& res_vocabs,
                      const RerankMode& mode,
                      const bool build_coarse_quantization,
                      const string& files_prefix,
                      const string& coarse_quantization_filename = "");
 private:
 /**
  * This function calculates for each point its main vocab and residual vocab
  * @param points_filename file with points in .fvecs or .bvecs format
  * @param points_count how many points should we handle
  * @param main_vocabs vocabularies for initial points
  * @param res_vocabs vocabularies for residuals
  */
  void PrepareCoarseQuantization(const string& points_filename,
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
  * @param transposed_coarse_quantizations result
  */
  void GetCoarseQuantizationsForSubset(const string& points_filename,
                                       const int start_pid,
                                       const int subset_size,
                                       const Centroids& main_vocabs,
                                       const Centroids& res_vocabs,
                                       vector<vector<ClusterId> >*
                                       transposed_coarse_quantizations);
 /**
  * This function serializes prepared coarse quantizations to file
  * @param transposed_coarse_quantizations quantizations to serialize.
  * They are transposed because of effective memory usage
  * @param filename file we should serialize to
  */
  void SerializeCoarseQuantizations(const vector<vector<ClusterId> >&
                                    transposed_coarse_quantizations,
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
  * This function reads point coarse quantization from file
  * @param pid identifier of target point
  * @param filename file with coarse quantizations
  * @param coarse_quantization result
  */
  void GetPointCoarseQuantization(const PointId pid,
                                  const string& filename,
                                  vector<ClusterId>* coarse_quantization);
 /**
  * This function restores counts of points from coarse quantizations
  * @param points_filename file with points in .fvecs or .bvecs format
  * @param points_count how many points should we index
  * @param main_vocabs vocabularies for initial points
  * @param res_vocabs vocabularies for residuals
  * We need them to init counts table correctly
  */
  void RestorePointsInCellsCountFromCourseQuantization(const string& points_filename,
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

 /**
  *  All index filenames will start from this prefix
  */
  string files_prefix_;
 /**
  *  Filename of file with coarse quantizations
  */
  string coarse_quantization_filename_;
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
};

template<class Record>
inline void GetRecord(const Point& point, const PointId pid,
                      const vector<ClusterId> coarse_quantization,
                      const vector<Centroids>& coarse_vocabs,
                      Record* result) {
}

template<class Record>
void InitParameters(const vector<Centroids>& fine_vocabs,
                    const RerankMode& mode,
                    const string& metainfo_filename) {
  gConfig.fine_vocabs = fine_vocabs;
  gConfig.rerank_mode = mode;
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
            Point temp(main_vocabs[0].size());
            cblas_saxpy(temp.size(), 1, &(main_vocabs[i][0]), 1, &(temp[0]), 1);
            cblas_saxpy(temp.size(), 1, &(res_vocabs[j][0]), 1, &(temp[0]), 1);
            precomputed_norms_[i][j] = cblas_sdot(temp.size(), &(temp[0]), 1, &(temp[0]), 1);
        }
    }
}


template<class Record>
void Indexer<Record>::SerializeCoarseQuantizations(const vector<vector<ClusterId> >&
		                                                          transposed_coarse_quantizations,
                                                        const string& filename) {
  ofstream quantizations_stream;
  quantizations_stream.open(filename.c_str(), ios::binary);
  if(!quantizations_stream.good()) {
    throw std::logic_error("Bad input stream");
  }
  cout << "Writing coarse quantizations started " << filename << endl;
  for(PointId pid = 0; pid < transposed_coarse_quantizations[0].size(); ++pid) {
    for(int index = 0; index < transposed_coarse_quantizations.size(); ++index) {
      ClusterId quantization = transposed_coarse_quantizations[index][pid];
      quantizations_stream.write((char*)&quantization, sizeof(quantization));
    }
  }
  quantizations_stream.close();
  cout << "Writing coarse quantizations started" << endl;
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
  cout << "Finish hierindex serializing....\n";
}

template<class Record>
void Indexer<Record>::GetCoarseQuantizationsForSubset(const string& points_filename,
                                                      const int start_pid,
                                                      const int subset_size,
                                                      const Centroids& main_vocabs,
                                                      const Centroids& res_vocabs,
                                                      vector<vector<ClusterId> >*
                                                            transposed_coarse_quantizations) {
  ifstream point_stream;
  point_stream.open(points_filename.c_str(), ios::binary);
  if(!point_stream.good()) {
    throw std::logic_error("Bad input points stream");
  }
  // we assume points are stored in .fvecs or .bvecs format
  point_stream.seekg(start_pid * (GetInputCoordSizeof() * SPACE_DIMENSION + sizeof(Dimensions)), ios::beg);
  vector<ClusterId> coarse_quantization(2);
  for(int point_number = 0; point_number < subset_size; ++point_number) {
    if(point_number % 1000 == 0) {
      cout << "Getting coarse quantization, point # " << start_pid + point_number << endl;
    }
    Point current_point;
    ReadPoint(point_stream, &current_point);
    ////// GETTING QUANTIZATIONS
    vector<float> main_products(main_vocabs.size());
    vector<float> res_products(res_vocabs.size());
    vector<pair<float, int> > temp(main_products.size());
    for (int i = 0; i < main_vocabs_.size(); ++i) {
      main_products[i] = cblas_sdot(point.size(), &(current_point[0]), 1, &(main_vocabs[i][0]), 1);
      temp[i] = std::make_pair(main_norms_[i] / 2 - main_products[i], i);
    }
    for (int i = 0; i < res_vocabs_.size(); ++i) {
      res_products[i] = cblas_sdot(current_point.size(), &(current_point[0]), 1, &(res_vocabs[i][0]), 1);
    }
    std::sort(temp.begin(), temp.end());
    ClusterId opt_main = 0;
    ClusterId opt_res = 0;
    float opt_distance = 99999999;
    for(int i = 0; i < 16; ++i) {
    for(int j = 0; j < res_vocabs_.size(); ++j) {
      int main_cluster = temp[i].second;
      float distance = precomputed_norms_[main_cluster][j] - 2 * main_products[main_cluster] - 2 * res_products[j];
      if(distance < opt_distance) {
        opt_main = main_cluster;
        opt_res = j;
        opt_distance = distance;
      }
    }
  }
    transposed_coarse_quantizations->at(0)[start_pid + point_number] = opt_main;
    coarse_quantization[0] = opt_main;
    transposed_coarse_quantizations->at(1)[start_pid + point_number] = opt_res;
    coarse_quantization[1] = opt_res;
    //////
    //ClusterId nearest = GetNearestClusterId(current_point, main_vocabs, 0, current_point.size() - 1);
    //Point residual;
    //GetResidual(current_point, main_vocabs[nearest], &residual);
    //ClusterId res_nearest = GetNearestClusterId(residual, res_vocabs, 0, current_point.size() - 1);
    //transposed_coarse_quantizations->at(0)[start_pid + point_number] = nearest;
    //coarse_quantization[0] = nearest;
    //transposed_coarse_quantizations->at(1)[start_pid + point_number] = res_nearest;
    //coarse_quantization[1] = res_nearest;
    int global_index = point_in_cells_count_.GetCellGlobalIndex(coarse_quantization);
    cell_counts_mutex_.lock();
    ++(point_in_cells_count_.table[global_index]);
    cell_counts_mutex_.unlock();
  }
}

template<class Record>
void Indexer<Record>::PrepareCoarseQuantization(const string& points_filename,
                                                const int points_count,
                                                const Centroids& main_vocabs,
                                                const Centroids& res_vocabs) {
  // we use transposed quantizations for efficient memory usage
  vector<vector<ClusterId> > transposed_coarse_quantizations; 
  transposed_coarse_quantizations.resize(2);
  transposed_coarse_quantizations[0].resize(points_count);
  transposed_coarse_quantizations[1].resize(points_count);
  vector<int> multiindex_table_dimensions;
  multiindex_table_dimensions.push_back(main_vocabs.size());
  multiindex_table_dimensions.push_back(res_vocabs.size());
  point_in_cells_count_.Resize(multiindex_table_dimensions);
  cout << "Memory for coarse quantizations allocated" << endl;
  boost::thread_group index_threads;
  int thread_points_count = points_count / THREADS_COUNT;
  for(int thread_id = 0; thread_id < THREADS_COUNT; ++thread_id) {
    PointId start_pid = thread_points_count * thread_id;
    index_threads.create_thread(boost::bind(&Indexer::GetCoarseQuantizationsForSubset,
                                            this, points_filename, start_pid, thread_points_count,
                                            boost::cref(main_vocabs), boost::cref(res_vocabs),
                                            &transposed_coarse_quantizations));
  }
  index_threads.join_all();
  if(coarse_quantization_filename_.empty()) {
    coarse_quantization_filename_ = files_prefix_ + "_coarse_quantizations.bin";
  }
  cout << "Coarse quantizations are calculated" << endl;
  SerializeCoarseQuantizations(transposed_coarse_quantizations, coarse_quantization_filename_);
  cout << "Coarse quantizations are serialized" << endl;
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
void Indexer<Record>::GetPointCoarseQuantization(const PointId pid,
                                                 const string& filename,
                                                 vector<ClusterId>* coarse_quantization) {
  ifstream coarse_quantization_stream;
  coarse_quantization_stream.open(filename.c_str(), ios::binary);
  if(!coarse_quantization_stream.good()) {
    throw std::logic_error("Bad input coarse quantizations stream");
  }
  coarse_quantization_stream.seekg((long long)pid * sizeof(ClusterId) * 2, ios::beg);
  coarse_quantization_stream.read((char*)&(coarse_quantization->at(0)),
                                  sizeof(coarse_quantization->at(0)));
  coarse_quantization_stream.read((char*)&(coarse_quantization->at(1)),
                                  sizeof(coarse_quantization->at(1)));
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
  vector<ClusterId> coarse_quantization(2);
  GetPointCoarseQuantization(start_pid + point_number,
                             coarse_quantization_filename_,
                             &coarse_quantization);
  cell_counts_mutex_.lock();
  int current_written_count = points_written_in_index->GetValue(coarse_quantization);
  int pid_multiindex = multiindex_.cell_edges.GetValue(coarse_quantization) + current_written_count;
  multiindex_.multiindex[pid_multiindex].pid = start_pid + point_number;
  // TODO - add rerank info
  //GetRecord<Record>(current_point, start_pid + point_number,
  //                  coarse_quantization, coarse_vocabs, &(multiindex_.multiindex[pid_multiindex]));
  points_written_in_index->SetValue(current_written_count + 1, coarse_quantization);
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
void Indexer<Record>::RestorePointsInCellsCountFromCourseQuantization(const string& points_filename,
                                                                      const int points_count,
                                                                      const Centroids& main_vocabs,
                                                                      const Centroids& res_vocabs) {
  vector<int> dimensions;
  dimensions.push_back(main_vocabs.size());
  dimensions.push_back(res_vocabs.size());
  point_in_cells_count_.Resize(dimensions);
  ifstream coarse_quantization_stream;
  coarse_quantization_stream.open(coarse_quantization_filename_.c_str(), ios::binary);
  if(!coarse_quantization_stream.good()) {
    throw std::logic_error("Bad input coarse quantizations stream");
  }
  CoarseQuantization quantization(2);
  for(PointId pid = 0; pid < points_count; ++pid) {
    if(pid % 100000 == 0) {
      cout << pid << endl;
    }
    coarse_quantization_stream.read((char*)&(quantization[0]), sizeof(ClusterId));
    coarse_quantization_stream.read((char*)&(quantization[1]), sizeof(ClusterId));
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
                                     const bool build_coarse_quantization,
                                     const string& files_prefix,
                                     const string& coarse_quantization_filename) {
  //InitParameters<Record>(fine_vocabs, mode, metainfo_filename);
  //InitBlasStructures(coarse_vocabs);
  files_prefix_ = files_prefix;
  coarse_quantization_filename_ = coarse_quantization_filename;
  PrecomputeEffectiveCentroidsNorms(main_vocabs, res_vocabs);
  cout << "Norms are precomputed" << endl;
  if(build_coarse_quantization) {
    PrepareCoarseQuantization(points_filename, points_count, main_vocabs, res_vocabs);
  } else {
  RestorePointsInCellsCountFromCourseQuantization(points_filename,
                                                  points_count,
                                                  main_vocabs,
                                                  res_vocabs);
  }
  FillHierIndex(points_filename, points_count, main_vocabs, res_vocabs, mode);
  cout << "Hierindex created" << endl;
  SerializeHierIndexFiles();
  cout << "Hierindex serialized" << endl;
}

//template<class Record>
//void MultiIndexer<Record>::InitBlasStructures(const vector<Centroids>& coarse_vocabs) {
//  coarse_vocabs_matrices_.resize(coarse_vocabs.size());
//  coarse_centroids_norms_.resize(coarse_vocabs.size(), vector<float>(coarse_vocabs[0].size()));
//  for(int coarse_id = 0; coarse_id < coarse_vocabs_matrices_.size(); ++coarse_id) {
//    coarse_vocabs_matrices_[coarse_id] = new float[coarse_vocabs[0].size() * coarse_vocabs[0][0].size()];
//    for(int i = 0; i < coarse_vocabs[0].size(); ++i) {
//      Coord norm = 0;
//      for(int j = 0; j < coarse_vocabs[0][0].size(); ++j) {
//        coarse_vocabs_matrices_[coarse_id][coarse_vocabs[0][0].size() * i + j] = coarse_vocabs[coarse_id][i][j];
//        norm += coarse_vocabs[coarse_id][i][j] * coarse_vocabs[coarse_id][i][j];
//      }
//      coarse_centroids_norms_[coarse_id][i] = norm;
//    }
//  }
//}

template<>
inline void GetRecord<PointId> (const Point& point, const PointId pid,
                                const vector<ClusterId> coarse_quantization,
                                const vector<Centroids>& coarse_vocabs,
                                PointId* result) {
  *result = pid;
}

inline void FillAdcInfo(const Point& point, const PointId pid,
                        const vector<Centroids>& fine_vocabs,
                        char* result) {
  int subvectors_count = fine_vocabs.size();
  int subvector_dim = point.size() / subvectors_count;
  for(int subvector_index = 0; subvector_index < subvectors_count; ++subvector_index) {
    Dimensions start_dim = subvector_index * subvector_dim;
    Dimensions final_dim = start_dim + subvector_dim;
    *((FineClusterId*)result) = (FineClusterId)GetNearestClusterId(point, fine_vocabs[subvector_index],
			                                                       start_dim, final_dim);
    result += sizeof(FineClusterId);
  }
}

template<>
inline void GetRecord<RerankADC8> (const Point& point, const PointId pid,
                                   const vector<ClusterId> coarse_quantization,
                                   const vector<Centroids>& coarse_vocabs,
                                   RerankADC8* result) {
  result->pid = pid;
  char* rerank_info_ptr = (char*)result + sizeof(pid);
  if(gConfig.rerank_mode == USE_RESIDUALS) {
    Point residual;
    GetResidual(point, coarse_quantization, coarse_vocabs, &residual);
    FillAdcInfo(residual, pid, gConfig.fine_vocabs, rerank_info_ptr);
  } else if (gConfig.rerank_mode == USE_INIT_POINTS) {
    FillAdcInfo(point, pid, gConfig.fine_vocabs, rerank_info_ptr);
  }
}

template<>
inline void GetRecord<RerankADC16> (const Point& point, const PointId pid,
                                    const vector<ClusterId> coarse_quantization,
                                    const vector<Centroids>& coarse_vocabs,
                                    RerankADC16* result) {
  result->pid = pid;
  char* rerank_info_ptr = (char*)result + sizeof(pid);
  if(gConfig.rerank_mode == USE_RESIDUALS) {
    Point residual;
    GetResidual(point, coarse_quantization, coarse_vocabs, &residual);
    FillAdcInfo(residual, pid, gConfig.fine_vocabs, rerank_info_ptr);
  } else if (gConfig.rerank_mode == USE_INIT_POINTS) {
    FillAdcInfo(point, pid, gConfig.fine_vocabs, rerank_info_ptr);
  }
}

#endif




