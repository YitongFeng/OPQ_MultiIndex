// OPQ_MI_DLL.cpp : Defines the exported functions for the DLL application.
//

#include "searcher.h"
#include "indexer.h"
#include "perfomance_util.h"

#include <iostream>
#include <time.h>
#include <thread>

using std::endl;
using std::ofstream;

Dimensions SPACE_DIMENSION;
string coarse_vocabs_file;    // File with vocabularies for multiindex structure
string fine_vocabs_file;	// File with vocabularies for reranking
RerankMode mode;	// Reranking approach, should be USE_RESIDUALS or USE_INIT_POINTS
string index_files_prefix;	// Common prefix of all multiindex files
string queries_file;	// File with queries .bvec or .fvec
PointType query_point_type;	//Type, should be BVEC or FVEC
string groundtruth_file;	//File with groundtruth
int queries_count;	//Number of queries to search
bool do_rerank;	  //Should we rerank?
int neighbours_count;	// Number of neighbours to look over
string report_file;		//File to write report in
int subspaces_centroids_count;	// Number of nearest centroids for each group of dimensions to handle

string imgFold;  // The img path prefix
string data_file;
string result_path;
int data_count;
int use_originaldata;
string db_path_file;
string query_path_file;
int k;
bool has_gt;
string data_name;
string root_path;
int fine_vocab_size;

vector<pair<string, Point>> queries;
unordered_map<string, vector<pair<string, int>>> groundtruth;
vector<string> db_paths;
MultiSearcher<RerankADC8, PointId> searcher;
Points dataset;

int SetOptions() {
	data_name = "fv31test";
	root_path = "I:\\Hashing\\Code\\output\\";
	imgFold = "I:\\trademark\\31\\";
	SPACE_DIMENSION = 4992;
	data_count = 103524;
	queries_count = 100;
	k = 100;
	fine_vocab_size = 8;
	neighbours_count = 150;	// 100, 500, 1000
	subspaces_centroids_count = 100;
	has_gt = false;

	//THREADS_COUNT = 8;
	//multiplicity = 2;
	string index_path = root_path + data_name + "\\";	// "I:\\Hashing\\Code\\OPQ_Mindex_fyt\\index\\fvImgDb\\"
	result_path = root_path + data_name + "\\result\\";
	queries_file = index_path + data_name + "_query_NP.fvecs";
	query_path_file = index_path + data_name + "_query_paths.txt";

	coarse_vocabs_file = index_path + data_name + "_coarse.dat";
	fine_vocabs_file = index_path + data_name + "_fine.dat";
	index_files_prefix = index_path;
	report_file = result_path + data_name + "_report_";
	groundtruth_file = root_path + "gt_744.txt";	//************need to prepare gt files! *********
	db_path_file = index_path + data_name + "_paths.txt";
	data_file = index_path + data_name + "_base_NP.fvecs";
	use_originaldata = 0;	// 计算距离的时候用原数据还是用opq code


	do_rerank = 0;
	mode = USE_RESIDUALS;
	query_point_type = FVEC;	  // FVEC, BVEC

	return 0;
}

bool init(){
	SetOptions();
	cout << "Options are set ...\n";

	//// Load dataset
	///*vector<pair<string, Point>> db(data_count);
	//loadDbOrQuery(data_file, db_path_file, db, query_point_type);*/
	//Points dataset;
	////ReadPoints<float, Coord>(data_file, &dataset, data_count);	

	// Load paths
	ifstream data_path_fid(db_path_file);
	db_paths.resize(data_count);
	for (int i = 0; i < data_count; i++){
		data_path_fid >> db_paths[i];
	}

	if (has_gt){
		// Load ground truth
		loadGt(groundtruth_file, imgFold, groundtruth);
		cout << "Groundtruth is read ...\n";
	}

	searcher.Init(index_files_prefix, coarse_vocabs_file,
		fine_vocabs_file, mode,
		subspaces_centroids_count,
		do_rerank, use_originaldata);

	return true;
}

//bool release(){
//}

//bool query(char* filepath, vector<pair<string, float>>& query_result, bool proposal_flag, vector<int> class_choose){
//
//	if (queries.size() == 0){
//		cerr << "Query size is 0!" << endl;
//		return false;
//	}
//	std::ofstream out(report_file.c_str(), ios::app);
//	vector<DistanceToPoint> result;
//
//	////////// 提取CNN特征
//
//	////////// 提取fv特征 得到：
//	//vector<Point> features;
//	Point feature;
//
//
//	float time = 0.0;
//	clock_t start = clock();
//
//	vector<DistanceToPoint> res_tmp;
//	searcher.GetNearestNeighbours(feature, neighbours_count * 2, &res_tmp, dataset);
//
//	// Remove the duplicates in res_one_img
//	removeDublicates(res_tmp, db_paths);
//	sort(res_tmp.begin(), res_tmp.end(), comp_with_dis);
//	if (res_tmp.size() >= k)
//		result = vector<DistanceToPoint>(res_tmp.begin(), res_tmp.begin() + k);
//	else
//		cerr << "Not enough search results!" << endl;
//
//	query_result.resize(k);
//	for (int i = 0; i < result.size(); i++){
//		string img_path = db_paths[result[i].second];
//		query_result[i] = make_pair(img_path, result[i].first);
//	}
//
//	clock_t end = clock();
//	time += (end - start);
//
//	// *************** Save Result ***************
//
//	if (has_gt)
//		showImages(filepath, result_path, result, groundtruth[splitFileName(filepath)], db_paths);
//	else
//		showImages(filepath, result_path, result, db_paths);
//
//
//	// **************** Cal Recall ***************
//	//float recall = 0.0;
//	//for (auto it = result.begin(); it != result.end(); ++it) //queries_count
//	//{
//	//	string imgname = it->first;
//	//	if (groundtruth.find(imgname) != groundtruth.end()){
//	//		auto gt = groundtruth[imgname];
//	//		auto rr = result[imgname];
//	//		recall += GetRecall(k, gt, rr, db_paths);
//	//	}
//	//}
//	//recall = recall / queries_count;*/
//	time = time * 1.0 / CLOCKS_PER_SEC / queries_count;
//
//	out.setf(ios::fixed);
//	cout << " " << time << " #N_" << neighbours_count << " " << endl;
//	out << " total time" << time * queries_count << " mean time: " << time
//		<< " #N_neighbors: " << neighbours_count << " " << endl;
//
//	return true;
//
//}

bool query(pair<string, Point> query_point, vector<pair<string, float>>& query_result, bool proposal_flag, int thread_num){

	if (queries.size() == 0){
		cerr << "Query size is 0!" << endl;
		return false;
	}
	char str_thread_num[10];
	sprintf(str_thread_num, "%d", thread_num);
	std::ofstream out(report_file + string(str_thread_num) + ".txt");
	out << "thread num is : " << thread_num << endl;
	if (&query_point)
		out << "query size is: " << query_point.second.size() << endl;
	else
		out << "query is empty!" << endl;

	vector<DistanceToPoint> result;

	////////// 提取CNN特征

	////////// 提取fv特征 得到：
	//vector<Point> features;
	Point feature;


	float time = 0.0;
	clock_t start = clock();

	vector<DistanceToPoint> res_tmp;
	searcher.GetNearestNeighbours(query_point.second, neighbours_count * 2, &res_tmp, dataset);

	// Remove the duplicates in res_one_img
	removeDublicates(res_tmp, db_paths);
	sort(res_tmp.begin(), res_tmp.end(), comp_with_dis);
	if (res_tmp.size() >= k){
		result = vector<DistanceToPoint>(res_tmp.begin(), res_tmp.begin() + k);
		out << "Success! Result num is: " << result.size() << endl;
	}
	else{
		out << "Result num is: " << res_tmp.size() << endl;
		cout << "Not enough search results!" << endl;
	}

	query_result.resize(k);
	for (int i = 0; i < result.size(); i++){
		string img_path = db_paths[result[i].second];
		query_result[i] = make_pair(img_path, result[i].first);
	}

	clock_t end = clock();
	time += (end - start);

	// *************** Save Result ***************

	if (has_gt)
		showImages(query_point.first, result_path, result, groundtruth[splitFileName(query_point.first)], db_paths);
	else
		showImages(query_point.first, result_path, result, db_paths);


	// **************** Cal Recall ***************
	//float recall = 0.0;
	//for (auto it = result.begin(); it != result.end(); ++it) //queries_count
	//{
	//	string imgname = it->first;
	//	if (groundtruth.find(imgname) != groundtruth.end()){
	//		auto gt = groundtruth[imgname];
	//		auto rr = result[imgname];
	//		recall += GetRecall(k, gt, rr, db_paths);
	//	}
	//}
	//recall = recall / queries_count;*/
	time = time * 1.0 / CLOCKS_PER_SEC / queries_count;

	/*out.setf(ios::fixed);
	cout << " " << time << " #N_" << neighbours_count << " " << endl;
	out << " total time" << time * queries_count << " mean time: " << time
		<< " #N_neighbors: " << neighbours_count << " " << endl;*/

	return true;
}

int main(){
	init();
	// Load queries
	queries.resize(queries_count);
	loadDbOrQuery(queries_file, query_path_file, queries, query_point_type);

	const int thread_num = 10;

	thread t[thread_num];
	int test_num = 10;
	for (int i = 0; i < thread_num; i++){
		vector<pair<string, float>> query_result;
		auto func = static_cast<bool(*)(pair<string, Point>, vector<pair<string, float>>&, bool, int)>(query);
		t[i] = thread(func, queries[i], query_result, false, i);
	}
	for (int i = 0; i < thread_num; i++)
		t[i].join();

	/*int test_num = 10;
	for (int i = 0; i < test_num; i++){
		vector<pair<string, float>> query_result;
		query(queries[i], query_result, false, 0);
	}*/

	return 0;
}