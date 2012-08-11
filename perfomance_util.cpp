// Copyright 2012 Yandex Artem Babenko
#include "perfomance_util.h"

extern string report_file;

PerfTester::PerfTester() {
	report_file_ = report_file;
	current_points_count = 0;
	handled_queries_count = 0;
	cells_traversed = 0;
	nearest_subcentroids_time = 0;
	cache_init_time = 0;
	merger_init_time = 0;
	full_traversal_time = 0;
	cell_coordinates_time = 0;
	cell_edges_time = 0;
	residual_time = 0;
	refining_time = 0;
	full_search_time = 0;

	for(int i = 0; i < 21; ++i) {
		list_length_thresholds_.push_back(std::pow(2.0, i));
	}
	current_threshold_index_ = 0;
	list_length_times_.resize(list_length_thresholds_.size(), 0.0);
}

void PerfTester::ResetQuerywiseStatistic() {
	current_threshold_index_ = 0;
	current_points_count = 0;
}

void PerfTester::NextNeighbour() {
	++current_points_count;
	if(current_points_count >= list_length_thresholds_[current_threshold_index_]) {
	  clock_t current_time = clock();
		list_length_times_[current_threshold_index_] += current_time - search_start;
		++current_threshold_index_;
	}
}

void PerfTester::DoReport(std::ofstream& out) {
  out << "Queries count: "
      << handled_queries_count << endl;
	out << "Average cells count: "
      << (double)cells_traversed / handled_queries_count << endl;
	out << "Average nearest subcentroids getting time: "
      << (double)nearest_subcentroids_time / handled_queries_count << endl;
	out << "Average cache init time: "
      << (double)cache_init_time / handled_queries_count << endl;
	out << "Average merger init time: "
      << (double)merger_init_time / handled_queries_count << endl;
	out << "Average full traversal time: "
      << (double)full_traversal_time / handled_queries_count << endl;
	out << "Average cells coordinates getting time: "
      << (double)cell_coordinates_time / handled_queries_count << endl;
	out << "Average cell edges getting time: "
      << (double)cell_edges_time/ handled_queries_count << endl;
	out << "Average residual time: "
      << (double)residual_time / handled_queries_count << endl;
	out << "Average refining time: "
      <<(double)refining_time / handled_queries_count << endl;
	out << "Average full search time: "
      << (double)full_search_time / handled_queries_count << endl;
}

void PerfTester::DoReport() {
	std::ofstream out(report_file_.c_str());
	DoReport(out);
}

int GetRecallAt(const int length, const vector<PointId>& groundtruth,
	              const vector<DistanceToPoint>& result) {
	if(groundtruth.empty()) {
    cout << "Groundtruth is empty!" << endl;
    return 0;
	}
	for(int index = 0; index < length && index < result.size(); ++index) {
		if(result[index].second == groundtruth[0]) {
      return 1;
		}
	}
	return 0;
}

double GetPresicionAt(const int length, const set<PointId>& groundtruth,
	                    const vector<DistanceToPoint>& result) {
	int found = 0;
	for(int index = 0; index < length && index < result.size() ; ++index) {
		if(groundtruth.find(result[index].second) != groundtruth.end()) {
		  found += 1;
		}
	}
	return (double)found / length; 
}