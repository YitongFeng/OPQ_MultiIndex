// Copyright 2012 Yandex Artem Babenko
#include <iostream>

//#include <boost/program_options.hpp>
#include <time.h>
#include <thread>
//#include <sys/time.h>
//#include <mkl.h>

#include "searcher.h"
#include "indexer.h"
#include "perfomance_util.h"


//#pragma comment(lib,"nearest_search_lib.lib")

//using namespace boost::program_options;
using std::endl;
using std::ofstream;

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

string imgFold;  // The img path prefix

string data_file;
string result_path;
int data_count;
int use_originaldata;
string db_path_file;
string query_path_file;
int k;
bool has_gt;

int SetOptions(int argc, char** argv) {
	string data_name = "CNN3000";
	string root_path = "I:\\Hashing\\Code\\output\\";
	imgFold = "T:\\整理的数据集\\";
	SPACE_DIMENSION = 256;
	data_count = 3745;
	queries_count = 100;
	k = 20;
	has_gt = false;

	//THREADS_COUNT = 8;
	//multiplicity = 2;
	string index_path = root_path + data_name + "\\";	// "I:\\Hashing\\Code\\OPQ_Mindex_fyt\\index\\fvImgDb\\"
	result_path = root_path + data_name + "\\results_thread\\";
	queries_file = root_path + "CNN_query_NP.fvecs";
	query_path_file = root_path + "CNN_query_paths.txt";

	coarse_vocabs_file = index_path + data_name + "_coarse.dat";
	fine_vocabs_file = index_path + data_name + "_fine.dat";
	index_files_prefix = index_path;
	report_file = result_path + data_name + "_report.txt";
	groundtruth_file = root_path + "gt_744.txt";	//************need to prepare gt files! *********
	db_path_file = index_path + data_name + "_paths.txt";
	

	neighbours_count = 100;	// 100, 500, 1000
	subspaces_centroids_count = 2 * k;
	data_file = index_path + data_name + "_base_NP.fvecs";
	use_originaldata = 0;	// 计算距离的时候用原数据还是用opq code


	do_rerank = 0;
	mode = USE_RESIDUALS;
	query_point_type = FVEC;	  // FVEC, BVEC

	return 0;
}


//template<class TSearcher>
//void TestSearcher(TSearcher& searcher,
//	const Points& queries,
//	const vector<vector<PointId> >& groundtruth, const Points& dataset) {
//	searcher.Init(index_files_prefix, coarse_vocabs_file,
//		fine_vocabs_file, mode,
//		subspaces_centroids_count,
//		do_rerank, use_originaldata); 
//
//	std::ofstream out(report_file.c_str(), ios::app);
//
//	vector< vector< DistanceToPoint> > result;
//	result.resize(queries_count);
//	
//	float time = 0;
//	clock_t start = clock();
//	for (int i = 0; i < queries_count; ++i) 
//	{
//		searcher.GetNearestNeighbours(queries[i], neighbours_count, &result[i], dataset);
//	}
//	clock_t end = clock();
//	time += (end - start) / CLOCKS_PER_SEC;
//	//gettimeofday(&end, NULL);
//	//time += diff_timeval(end, start);
//
//	// **************** Cal Recall ***************
//	float recall = 0.0;
//	for (int i = 0; i < queries_count; ++i) //queries_count
//	{
//		recall += GetRecall(k, groundtruth[i], result[i]);
//	}
//
//	recall = recall / queries_count;
//	time = time / queries_count;
//
//	out.setf(ios::fixed);
//	cout << recall << " " << time << " #N_" << neighbours_count << " " << endl;
//	out << recall << " " << time << " #N_" << neighbours_count << " " << endl;
//}

template<class TSearcher>
void TestSearcher(TSearcher& searcher,
	vector<pair<string, Point>>& queries,
	unordered_map<string, vector<pair<string, int>>> groundtruth, Points& dataset, vector<string>& db_paths) {
	
	searcher.Init(index_files_prefix, coarse_vocabs_file,
		fine_vocabs_file, mode,
		subspaces_centroids_count,
		do_rerank, use_originaldata); 

	std::ofstream out(report_file.c_str(), ios::app);

	if (queries.size() == 0){
		cerr << "Query size is 0!" << endl;
		return;
	}
	unordered_map<ImageName, vector< DistanceToPoint> > result;

	string cur_img = queries[0].first;
	float time = 0.0;
	clock_t start = clock();
	vector<DistanceToPoint> res_one_img;

	for (int i = 0; i < queries_count; ++i) 
	{
		vector<DistanceToPoint> res_tmp;
		searcher.GetNearestNeighbours(queries[i].second, neighbours_count, &res_tmp, dataset);

		if (queries[i].first != cur_img){
			// Remove the duplicates in res_one_img
			removeDublicates(res_one_img, db_paths);
			sort(res_one_img.begin(), res_one_img.end(), comp_with_dis);
			string cur_name = splitFileName(cur_img);
			if (res_one_img.size() >= k)
				result[cur_name] = vector<DistanceToPoint>(res_one_img.begin(), res_one_img.begin() + k);
			else
				cerr << "Not enough retrieval results! " << endl;
			res_one_img.clear();
			cur_img = queries[i].first;
		}
		res_one_img.insert(res_one_img.end(), res_tmp.begin(), res_tmp.end());
	}
	if (!res_one_img.empty()){
		// Remove the duplicates in res_one_img
		removeDublicates(res_one_img, db_paths);
		sort(res_one_img.begin(), res_one_img.end(), comp_with_dis);
		string cur_name = splitFileName(cur_img);
		if (res_one_img.size() >= neighbours_count)
			result[cur_name] = vector<DistanceToPoint>(res_one_img.begin(), res_one_img.begin() + neighbours_count);
		else
			cerr << "Not enough retrieval results! " << endl;
		res_one_img.clear();
	}
	clock_t end = clock();
	time += (end - start);

	// *************** Save Result ***************
	for (auto it = result.begin(); it != result.end(); it++){
		ImagePath query_path = imgFold + "Img_Db\\" + it->first;
		if (has_gt)
			showImages(query_path, result_path, it->second, groundtruth[it->first], db_paths);
		else
			showImages(query_path, result_path, it->second, db_paths);
	}

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

	/*recall = recall / queries_count;*/
	time = time * 1.0 / CLOCKS_PER_SEC/ queries_count;

	out.setf(ios::fixed);
	cout << " " << time << " #N_" << neighbours_count << " " << endl;
	out << " total time" << time * queries_count  << " mean time: " << time 
		<< " #N_neighbors: " << neighbours_count << " " << endl;
}


int main(int argc, char** argv) {
	SetOptions(argc, argv);
	cout << "Options are set ...\n";

	// Load queries
	vector<pair<string, Point>> queries(queries_count);
	loadDbOrQuery(queries_file, query_path_file, queries, query_point_type);

	// Load dataset
	/*vector<pair<string, Point>> db(data_count);
	loadDbOrQuery(data_file, db_path_file, db, query_point_type);*/
	Points dataset;
	//ReadPoints<float, Coord>(data_file, &dataset, data_count);	

	// Load paths
	ifstream data_path_fid(db_path_file);
	vector<string> db_paths(data_count);
	for (int i = 0; i < data_count; i++){
		data_path_fid >> db_paths[i];
	}

	unordered_map<string, vector<pair<string, int>>> groundtruth;
	if (has_gt){
		// Load ground truth
		loadGt(groundtruth_file, imgFold, groundtruth);
		cout << "Groundtruth is read ...\n";
	}

	vector<Centroids> fine_vocabs;
	ReadFineVocabs<float>(fine_vocabs_file, &fine_vocabs);

	const int thread_num = 2;
	thread t[thread_num];

	if (fine_vocabs.size() == 8) {
		MultiSearcher<RerankADC8, PointId> searcher;
		for (int i = 0; i < thread_num; i++){
			vector<pair<string, Point>> tmp_query(1, queries[i]);
			auto func = static_cast<void(*)(MultiSearcher<RerankADC8, PointId>&, vector<pair<string, Point>>&, unordered_map<string, vector<pair<string, int>>>,
				Points&, vector<string>&) > (TestSearcher);
			t[i] = thread(func, ref(searcher), ref(tmp_query), groundtruth, ref(dataset), db_paths);

		}
		for (int i = 0; i < thread_num; i++){
			t[i].join();
		}
		//TestSearcher<MultiSearcher<RerankADC8, PointId> >(searcher, queries, groundtruth, dataset, db_paths);
		searcher.GetPerfTester().DoReport();
		//TestSearcher<MultiSearcher<RerankADC8, PointId> >(searcher, queries, groundtruth, dataset, paths);
	}
	else if (fine_vocabs.size() == 16) {
		MultiSearcher<RerankADC16, PointId> searcher;
		TestSearcher<MultiSearcher<RerankADC16, PointId> >(searcher, queries, groundtruth, dataset, db_paths);
		searcher.GetPerfTester().DoReport();
		//TestSearcher<MultiSearcher<RerankADC16, PointId> >(searcher, queries, groundtruth, dataset, paths);
	}

	return 0;
}
