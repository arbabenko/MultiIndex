// Copyright 2012 Yandex Artem Babenko
#include <iostream>

#include <boost/program_options.hpp>

#include "indexer.h"

using namespace boost::program_options;

/**
 * Number of threads for indexing
 */
int THREADS_COUNT;
/**
 * Type, should be BVEC or FVEC
 */
PointType point_type;
/**
 * Number of coordinates in a point
 */
Dimensions SPACE_DIMENSION;
/**
 * File with vocabularies for main centroids
 */
string main_vocabs_file;
/**
 * File with vocabularies for residuals
 */
string res_vocabs_file;
/**
 * File with points to index
 */
string points_file;
/**
 * File with points metainfo (imageId, etc.)
 */
string metainfo_file;
/**
 * Reranking approach, should be USE_RESIDUALS or USE_INIT_POINTS
 */
RerankMode mode;
/**
 * Common prefix of all multiindex files
 */
string files_prefix;
/**
 * Should we calculate coarse quantizations (they can be precomputed)
 */
bool build_quantizations;
/**
 * File with points coarse quantizations
 */
string quantizations_file;
/**
 * How many points should we index
 */
int points_count;
/**
 * Number of main centroids to handle for indexing
 */
int main_centroids_count;

int SetOptions(int argc, char** argv) {
  options_description description("Options");
  description.add_options()
    ("threads_count,t", value<int>())
    ("points_file,p", value<string>())
    ("metainfo_file,z", value<string>())
    ("main_vocabs_file,c", value<string>())
    ("res_vocabs_file,f", value<string>())
    ("input_point_type,i", value<string>())
    ("build_quantizations,b", bool_switch(), "Flag B")
    ("use_residuals,r", bool_switch(), "Flag R")
    ("points_count,p", value<int>())
    ("quantization_file,q", value<string>())
    ("space_dim,d", value<int>())
    ("main_centroids_count,s", value<int>())
    ("files_prefix,_", value<string>());
  variables_map name_to_value;
  try {
    store(command_line_parser(argc, argv).options(description).run(), name_to_value);
  } catch (const invalid_command_line_syntax& inv_syntax) {
    switch (inv_syntax.kind()) {
      case invalid_syntax::missing_parameter :
        cout << "Missing argument for option '" << inv_syntax.tokens() << "'.\n";
        break;
      default:
        cout << "Syntax error, kind " << int(inv_syntax.kind()) << "\n";
        break;
      };
    return 1;
  } catch (const unknown_option& unkn_option) {
    cout << "Unknown option '" << unkn_option.get_option_name() << "'\n";
    return 1;
  }
  if (name_to_value.count("help")) {
    cout << description << "\n";
    return 1;
  }

  THREADS_COUNT =              name_to_value["threads_count"].as<int>();
  points_file =                name_to_value["points_file"].as<string>();
  metainfo_file =              name_to_value["metainfo_file"].as<string>();
  main_vocabs_file =           name_to_value["main_vocabs_file"].as<string>();
  res_vocabs_file =            name_to_value["res_vocabs_file"].as<string>();
  SPACE_DIMENSION =            name_to_value["space_dim"].as<int>();
  files_prefix =               name_to_value["files_prefix"].as<string>();
  points_count =               name_to_value["points_count"].as<int>();
  main_centroids_count =       name_to_value["main_centroids_count"].as<int>();
 
  build_quantizations = (name_to_value["build_quantizations"].as<bool>() == true) ? true : false;
  mode = name_to_value["use_residuals"].as<bool>() == true ? USE_RESIDUALS : USE_INIT_POINTS;

  if (name_to_value.find("quantization_file") != name_to_value.end()) {
    quantizations_file =  name_to_value["quantization_file"].as<string>();
  }
  if (name_to_value["input_point_type"].as<string>() == "FVEC") {
    point_type = FVEC;
  } else if(name_to_value["input_point_type"].as<string>() == "BVEC") {
    point_type = BVEC;
  }
  return 0;
}

int main(int argc, char** argv) {
  SetOptions(argc, argv);
  cout << "Options are set ...\n";
  Centroids main_vocabs;
  Centroids res_vocabs;
  ReadCentroids<float>(main_vocabs_file, SPACE_DIMENSION, &main_vocabs);
  ReadCentroids<float>(res_vocabs_file, SPACE_DIMENSION, &res_vocabs);
  cout << "Vocs are read ...\n";
  Indexer<RerankADC8> indexer;
  indexer.BuildHierIndex(points_file, metainfo_file, points_count,
                         main_vocabs, res_vocabs, mode,
                         build_quantizations, files_prefix,
                         quantizations_file, main_centroids_count);
  return 0;
}