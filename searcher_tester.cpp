// Copyright 2012 Yandex Artem Babenko
#include <iostream>

#include <boost/program_options.hpp>

#include <mkl.h>

#include "searcher.h"
#include "indexer.h"

using namespace boost::program_options;

/**
 * Number of threads for indexing
 */
Dimensions SPACE_DIMENSION;
/**
 * File with vocabularies for multiindex structure
 */
string coarse_vocabs_file;
/**
 * File with vocabularies for reranking
 */
string fine_vocabs_file;
/**
 * Reranking approach, should be USE_RESIDUALS or USE_INIT_POINTS
 */
RerankMode mode;
/**
 * Common prefix of all multiindex files
 */
string index_files_prefix;
/**
 * File with queries (.bvec or .fvec)
 */
string queries_file;
/**
 * Type, should be BVEC or FVEC
 */
PointType query_point_type;
/**
 * File with groundtruth (.ivec)
 */
string groundtruth_file;
/**
 * Number of queries to search
 */
int queries_count;
/**
 * Should we rerank?
 */
bool do_rerank;
/**
 * Number of neighbours to look over
 */
int neighbours_count;
/**
 * File to write report in
 */
string report_file;
/**
 * Number of nearest centroids for each group of dimensions to handle
 */
int subspaces_centroids_count;



int SetOptions(int argc, char** argv) {
  options_description description("Options");
  description.add_options()
    ("index_files_prefix,i", value<string>())
    ("queries_file,q", value<string>())
    ("queries_count,n", value<int>())
    ("neighbours_count,k", value<int>())
    ("groundtruth_file,g", value<string>())
    ("coarse_vocabs_file,c", value<string>())
    ("fine_vocabs_file,f", value<string>())
    ("query_point_type,t", value<string>())
    ("do_rerank,l", bool_switch(), "Flag B")
    ("use_residuals,r", bool_switch(), "Flag R")
    ("points_count,p", value<int>())
    ("report_file,o", value<string>())
    ("space_dim,d", value<int>())
    ("subspaces_centroids_count,s", value<int>());
  variables_map name_to_value;
  try {
    store(command_line_parser(argc, argv).options(description).run(), name_to_value);
  } catch (const invalid_command_line_syntax &inv_syntax) {
    switch (inv_syntax.kind()) {
      case invalid_syntax::missing_parameter :
        cout << "Missing argument for option '" << inv_syntax.tokens() << "'.\n";
        break;
      default:
        cout << "Syntax error, kind " << int(inv_syntax.kind()) << "\n";
        break;
       };
    return 1;
  } catch (const unknown_option &unkn_opt) {
    cout << "Unknown option '" << unkn_opt.get_option_name() << "'\n";
    return 1;
  }
  if (name_to_value.count("help")) {
    cout << description << "\n";
    return 1;
  }

  coarse_vocabs_file =         name_to_value["coarse_vocabs_file"].as<string>();
  fine_vocabs_file =           name_to_value["fine_vocabs_file"].as<string>();
  SPACE_DIMENSION =            name_to_value["space_dim"].as<int>();
  index_files_prefix =         name_to_value["index_files_prefix"].as<string>();
  queries_file =               name_to_value["queries_file"].as<string>();
  report_file =                name_to_value["report_file"].as<string>();
  groundtruth_file =           name_to_value["groundtruth_file"].as<string>();
  queries_count =              name_to_value["queries_count"].as<int>();
  neighbours_count =           name_to_value["neighbours_count"].as<int>();
  subspaces_centroids_count =  name_to_value["subspaces_centroids_count"].as<int>();
 
  do_rerank =                  (name_to_value["do_rerank"].as<bool>() == true) ? true : false;
  mode =                       (name_to_value["use_residuals"].as<bool>() == true) ? USE_RESIDUALS : USE_INIT_POINTS;
  if (name_to_value["query_point_type"].as<string>() == "FVEC") {
    query_point_type = FVEC;
  } else if(name_to_value["query_point_type"].as<string>() == "BVEC") {
    query_point_type = BVEC;
  }
  return 0;
}

template<class TSearcher>
void TestSearcher(TSearcher& searcher,
                  const Points& queries,
                  const vector<vector<PointId> >& groundtruth) {
  searcher.Init(index_files_prefix, coarse_vocabs_file,
                fine_vocabs_file, mode,
                subspaces_centroids_count,
                do_rerank);
  cout << "Searcher inited ...\n";
  vector<DistanceToPoint> result;
  vector<int> counts;
  for(int i = 0; i < 21; ++i){
      counts.push_back(std::pow(2.0,i));
  }
  //counts.push_back(neighbours_count);
  for(int k = 0; k < counts.size(); ++k) {
    float recall = 0.0;
    vector<double> recalls(5, 0.0);
    clock_t start = clock();
    for(int i = 0; i < queries_count; ++i) {
      std::cout << i << std::endl;
      neighbours_count = counts[k];
      result.clear();
      searcher.GetNearestNeighbours(queries[i], neighbours_count, &result);
      recalls[0] += GetRecallAt(1, groundtruth[i], result);
      recalls[1] += GetRecallAt(10, groundtruth[i], result);
      recalls[2] += GetRecallAt(100, groundtruth[i], result);
      recalls[3] += GetRecallAt(1000, groundtruth[i], result);
      recalls[4] += GetRecallAt(10000, groundtruth[i], result);
      recall += GetRecallAt(result.size(), groundtruth[i], result);
      result.clear();
    }
    //cout << "R@1 "     << recalls[0] / queries_count << "\n" <<
    //        "R@10 "    << recalls[1] / queries_count << "\n" <<
    //        "R@100 "   << recalls[2] / queries_count << "\n" <<
    //        "R@1000 "  << recalls[3] / queries_count << "\n" <<
    //        "R@10000 " << recalls[4] / queries_count << endl;
    cout << "R@" << neighbours_count << " " << recall / queries_count << endl;
    searcher.GetPerfTester().DoReport();
    clock_t finish = clock();
    //std::cout << "Average search time(ms): "<<(double)(finish - start) / queries.size() << std::endl;
  }
}

int main(int argc, char** argv) {
  SetOptions(argc, argv);
  cout << "Options are set ...\n";
  Points queries;
  if(query_point_type == BVEC) {
    ReadPoints<unsigned char, Coord>(queries_file, &queries, queries_count);
  } else if (query_point_type == FVEC) {
    ReadPoints<float, Coord>(queries_file, &queries, queries_count);
  }
  cout << "Queries are read ...\n";
  vector<vector<PointId> > groundtruth;
  ReadPoints<int, PointId>(groundtruth_file, &groundtruth, queries_count);
  MKL_Set_Num_Threads(1);
  cout << "Groundtruth is read ...\n";
  vector<Centroids> fine_vocabs;
  ReadFineVocabs<float>(fine_vocabs_file, &fine_vocabs);
  if(fine_vocabs.size() == 8) {
    MultiSearcher<RerankADC8, PointId> searcher;
    TestSearcher<MultiSearcher<RerankADC8, PointId> > (searcher, queries, groundtruth);
  } else if(fine_vocabs.size() == 16) {
    MultiSearcher<RerankADC16, PointId> searcher;
    TestSearcher<MultiSearcher<RerankADC16, PointId> > (searcher, queries, groundtruth);
  }
  return 0;
}
