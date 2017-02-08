/** @file */
// Copyright 2012 Yandex Artem Babenko
#ifndef _PERFOMANCE_UTIL_H_
#define _PERFOMANCE_UTIL_H_

#include <iostream>
#include <set>
#include <vector>
#include <io.h>
#include <opencv2\opencv.hpp>
#include <opencv2\core\core.hpp>

#include "data_util.h"


using std::cout;
using std::endl;
using std::ofstream;
using std::pair;
using std::set;
using std::vector;
using cv::Mat;

extern string report_file;

/**
* \typedef
*  Typedef for point identifier and distance from query
*/
typedef pair<Distance, PointId> DistanceToPoint;



/**
*  This simple class stores timing of search working process
*/
class PerfTester {
public:
	PerfTester();
	/**
	*  Number of neighbours already found
	*/
	int current_points_count;
	/**
	*  Pretty report of timing
	*/
	void DoReport();
	/**
	*  Reset all prevoius statistic before
	*  new query handling
	*/
	void ResetQuerywiseStatistic();
	/**
	*  Signal about next point
	*/
	void NextNeighbour();
	/**
	*  Number of handled queries
	*/
	int handled_queries_count;
	/**
	*  Number of traversed items of multiindex
	*/
	int cells_traversed;
	unsigned long long nearest_subcentroids_time;
	unsigned long long cache_init_time;
	unsigned long long merger_init_time;
	unsigned long long full_traversal_time;
	unsigned long long cell_coordinates_time;
	unsigned long long cell_edges_time;
	unsigned long long residual_time;
	unsigned long long refining_time;
	unsigned long long full_search_time;
	unsigned long long search_start;
private:
	string report_file_;
	void DoReport(ofstream& out);
	vector<int> list_length_thresholds_;
	int current_threshold_index_;
	vector<float> list_length_times_;
};

inline bool comp_with_dis(DistanceToPoint r1, DistanceToPoint r2){
	return r1.first < r2.first;
};

inline bool comp_with_pointid(DistanceToPoint r1, DistanceToPoint r2){
	return r1.second < r2.second;
};

// From the result_tmp, remove the result with the same image path
void removeDublicates(vector<DistanceToPoint>& res_tmp, vector<string>& db_paths);

// Save images of searching results
void showImages(ImagePath queryName, string saveFold, vector< DistanceToPoint>& qRe, vector<pair<string, int> > qOrder, vector<string> paths);
void showImages(ImagePath queryName, string saveFold, vector< DistanceToPoint>& qRe, vector<string>& paths);
void showImages(ImagePath queryName, string saveFold, vector< DistanceToPoint>& qRe, vector<string>& cnn_paths, vector<string>& fv_paths);

/**
*  This function returns recall at specified length
* @param length specified size of search results
* @param groundtruth groundtruth
* @param result search results
*/
int GetRecallAt(const int length, const vector<PointId>& groundtruth,
	const vector<DistanceToPoint>& result);
/**
*  This function returns precision at specified length
* @param length specified size of search results
* @param groundtruth groundtruth
* @param result search results
*/
double GetPresicionAt(const int length, const vector<PointId>& groundtruth,
	const vector<DistanceToPoint>& result);

/**
*  This function returns recall at full length
*/
double GetRecall(int k, const vector<PointId>& groundtruth,
	const vector<DistanceToPoint>& result);
double GetRecall(int k, vector<pair<string, int>>& gt_i,
	vector<DistanceToPoint>& res_i, vector<string>& db_paths);

float get_distance(Point data, Point query);

float compute_relative_distance_error(int k, const Points& dataset, const Point& query, const vector<PointId>& groundtruth, const vector<DistanceToPoint>& result);

float compute_number_closer(int k, const vector<PointId>& groundtruth, const vector<DistanceToPoint>& result);

float compute_mean_reciprocal_rank(int k, const vector<PointId>& groundtruth, const vector<DistanceToPoint>& result);

float compute_mean_average_precision(int k, const vector<PointId>& groundtruth, const vector<DistanceToPoint>& result);

float compute_discounted_culmulative_gain(int k, const vector<PointId>& groundtruth, const vector<DistanceToPoint>& result);


#endif
