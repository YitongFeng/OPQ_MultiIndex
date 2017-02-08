// OPQ_MI_DLL.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include "OPQ_MI_DLL.h"
#include "searcher.h"
#include "indexer.h"
#include "perfomance_util.h"
#include "sift_fv_dll.h"
#include "preprocess_dll.h"
#include "CNN_API.hpp"

#include <iostream>
#include <time.h>
#include <thread>

using std::endl;
using std::ofstream;

Dimensions FV_DIMENSION;
Dimensions CNN_DIMENSION;	 
string coarse_vocabs_file_fv;    // File with vocabularies for multiindex structure
string coarse_vocabs_file_cnn;
string fine_vocabs_file_fv;	// File with vocabularies for reranking
string fine_vocabs_file_cnn;
RerankMode mode;	// Reranking approach, should be USE_RESIDUALS or USE_INIT_POINTS
string index_files_prefix;	// Common prefix of all multiindex files
string queries_file;	// File with queries (.bvec or .fvec)
PointType query_point_type;	//Type, should be BVEC or FVEC
string groundtruth_file;	//File with groundtruth
int queries_count;	//Number of queries to search
bool do_rerank;	  //Should we rerank?
int neighbours_count;	// Number of neighbours to look over
string report_file;		//File to write report in
int subspaces_centroids_count;	// Number of nearest centroids for each group of dimensions to handle
int cnn_res_num;	// how many results of cnn keep

string imgFold;  // The img path prefix
string data_file;
string result_path;
int data_count_fv;
int data_count_cnn;
int use_originaldata;
string db_path_file_fv;
string db_path_file_cnn;
string query_path_file;
int k;
bool has_gt;
string data_name_cnn;
string data_name_fv;
string root_path;
string pcaM_path;
string cnnR_path;
string fvR_path;
int fine_vocab_size;

vector<pair<string, Point>> queries;
unordered_map<string, vector<pair<string, int>>> groundtruth;
vector<string> db_paths_fv;
vector<string> db_paths_cnn;
MultiSearcher<RerankADC8, PointId> searcher;
Points dataset;
Mat pcaM;
Mat cnnR;
Mat fvR;

//#define debug
//#define debug_fv	// 开启的时候要注释掉关于cnn的一些代码
//#define debug_cnn	// 开启的时候要注释掉一些fv的代码

#ifdef debug
	ofstream f_debug("debug.txt");
#endif

int SetOptions() {
	data_name_fv = "fv_final2";
	data_name_cnn = "CNN2w";
	root_path = "Data\\";
	imgFold = "I:/trademark/31/";	//只有load gt的时候用
	data_count_fv = 227430;
	data_count_cnn = 225533;
	
	CNN_DIMENSION = 256;
	FV_DIMENSION = 1024;	
	queries_count = 100;
	k = 100;
	fine_vocab_size = 8;
	neighbours_count = 150;	// 100, 500, 1000
	subspaces_centroids_count = 40;
	cnn_res_num = 10;
	has_gt = false;

	//THREADS_COUNT = 8;
	//multiplicity = 2;
	string index_path = root_path;	// "I:\\Hashing\\Code\\OPQ_Mindex_fyt\\index\\fvImgDb\\"
	result_path = root_path + "result\\";
	//queries_file = root_path + "CNN_query_NP.fvecs";
	//query_path_file = root_path + "CNN_query_paths.txt";
	pcaM_path = index_path + data_name_fv + "_pca1024.xml";	
	cnnR_path = index_path + data_name_cnn + "_R.xml";
	fvR_path = index_path + data_name_fv + "_R.xml";
	coarse_vocabs_file_fv = index_path + data_name_fv + "_coarse.dat";
	fine_vocabs_file_fv = index_path + data_name_fv + "_fine.dat";
	coarse_vocabs_file_cnn = index_path + data_name_cnn + "_coarse.dat";
	fine_vocabs_file_cnn = index_path + data_name_cnn + "_fine.dat";
	index_files_prefix = index_path;
	report_file = result_path + "_report.txt";
	groundtruth_file = root_path + "gt_744.txt";	//************need to prepare gt files! *********
	db_path_file_fv = index_path + data_name_fv + "_paths.txt";
	db_path_file_cnn = index_path + data_name_cnn + "_paths.txt";
	
	use_originaldata = 0;	// 计算距离的时候用原数据还是用opq code


	do_rerank = 0;
	mode = USE_RESIDUALS;
	query_point_type = FVEC;	  // FVEC, BVEC

	return 0;
}

OPQ_MI_DLL bool init(){
	SetOptions();
	cout << "Options are set ...\n";

	init_preprocess();
#ifndef debug_cnn
	init_fv();
#endif
#ifndef debug_fv
	cnn_init();
#endif
	// ********* Load pca ***************
	cv::FileStorage fs(pcaM_path.c_str(), cv::FileStorage::READ);
	fs["eigenV"] >> pcaM;
	fs.release();
	// ********** Load R ****************
	cv::FileStorage fscr(cnnR_path.c_str(), cv::FileStorage::READ);
	fscr["eigenV"] >> cnnR;
	fscr.release();

	cv::FileStorage fsvr(fvR_path.c_str(), cv::FileStorage::READ);
	fsvr["eigenV"] >> fvR;
	fsvr.release();

#ifdef debug
	for (int i = 0; i < 10; i++){
		f_debug << fvR.at<float>(0, i) << " ";
	}
	f_debug << endl << "above is fv R" << endl;
#endif 

	// Load paths
	ifstream data_path_fid_fv(db_path_file_fv);
	db_paths_fv.resize(data_count_fv);
	for (int i = 0; i < data_count_fv; i++){
		data_path_fid_fv >> db_paths_fv[i];
	}

	ifstream data_path_fid_cnn(db_path_file_cnn);
	db_paths_cnn.resize(data_count_cnn);
	for (int i = 0; i < data_count_cnn; i++){
		data_path_fid_cnn >> db_paths_cnn[i];
	}

	if (has_gt){
		// Load ground truth
		loadGt(groundtruth_file, imgFold, groundtruth);
		cout << "Groundtruth is read ...\n";
	}
	
	searcher.Init(index_files_prefix, coarse_vocabs_file_fv, coarse_vocabs_file_cnn,
		fine_vocabs_file_fv, fine_vocabs_file_cnn, mode, subspaces_centroids_count,
		do_rerank, use_originaldata);

#ifdef debug
	searcher.setDebugTrue();
#endif // debug

	return true;
}

//OPQ_MI_DLL bool release(){
//}

OPQ_MI_DLL bool query(const char* filepath, vector<pair<string, float>>& query_result){
	cout << "begin query..." << endl;
	
	Mat query_img = cv::imread(filepath);
	
	/*char str_thread_num[10];
	sprintf(str_thread_num, "%d", thread_num);
	std::ofstream out(report_file + string(str_thread_num) + ".txt");
	out << "thread num is : " << thread_num << endl;
	if (&query_img)
		out << "Load query image success!"  << endl;
	else
		out << "query is empty!" << endl;*/

	vector<DistanceToPoint> result;

	Mat image_edge_out;
	Mat image_out;
	preprocess(query_img, image_edge_out, image_out);
	cout << image_out.rows << endl;

	////////// 提取CNN特征
	cout << "extract cnn feature..." << endl;
	Mat cnn;
#ifndef debug_fv
	cnn_fea(image_out, cnn);
	cnn = cnn * cnnR.t();
#endif

	cout << "cnn fea sccess! fea size is: " << cnn.size() << endl;
	Point cnn_feature(CNN_DIMENSION, 0);

#ifndef debug_fv
	for (int i = 0; i < cnn.cols; i++){
		cnn_feature[i] = cnn.at<float>(0, i);
	}
#endif

	////////// 提取fv特征
	cout << " 提取fv特征..." << endl;
#ifndef debug_cnn
	Mat fv;
	compute_fv(image_out, image_edge_out, fv);
	cout << fv.cols << endl;

#ifdef debug
	f_debug << endl << "original fv feature: " << endl;
	for (int i = 0; i < 10; i++){
		f_debug << fv.at<float>(0, i) << " ";
	}
	/*f_debug << "\n after hog feature from i=511: " << endl;
	for (int i = 511; i < 522; i++)
		f_debug << fv.at<float>(0, i) << " ";*/
#endif

	cout << "Do pca..." << endl;
	fv = fv * pcaM.t();		// ******************* 需要修改 ************
	fv = fv * fvR.t();
	cout << "pca done!" << endl;

#endif
	Point fv_feature(FV_DIMENSION, 0);


#ifndef debug_cnn
	cout << "final fv fea: " << endl;
	for (int i = 0; i < fv.cols; i++){
		fv_feature[i] = fv.at<float>(0, i);
	}
#endif

#ifdef debug
	f_debug << "final fv fea: " << endl;
	for (int i = 0; i < 10; i++){
		if (i < 10)
			f_debug << fv_feature[i] << " ";
	}
#endif
	cout << "extract feature success! " << endl;
	pair<string, Point> query_point_fv = std::make_pair(string(filepath), fv_feature);
	pair<string, Point> query_point_cnn = std::make_pair(string(filepath), cnn_feature);


	float time = 0.0;
	clock_t start = clock();

	vector<DistanceToPoint> res_cnn, res_fv, res_tmp;
#ifndef debug_fv
	searcher.GetNearestNeighbours(query_point_cnn.second, neighbours_count, &res_cnn, dataset, 0);
	cout << "get cnn result success! " << endl;
	removeDublicates(res_cnn, db_paths_cnn);
	sort(res_cnn.begin(), res_cnn.end(), comp_with_dis);
	if (res_cnn.size() > cnn_res_num)
		res_tmp.insert(res_tmp.end(), res_cnn.begin(), res_cnn.begin() + cnn_res_num);
	else{
#ifdef debug
		f_debug << "not enough cnn result!" << endl;
#endif
		res_tmp.insert(res_tmp.end(), res_cnn.begin(), res_cnn.end());
	}
#endif
	searcher.GetNearestNeighbours(query_point_fv.second, neighbours_count * 3, &res_fv, dataset, 1);
	cout << "Get fv result success! " << endl;
	removeDublicates(res_fv, db_paths_fv);
	sort(res_fv.begin(), res_fv.end(), comp_with_dis);
	res_tmp.insert(res_tmp.end(), res_fv.begin(), res_fv.end());

	// TODO: 如果cnn和fv的db库一样，要修改。

	if (res_tmp.size() >= k){
		result = vector<DistanceToPoint>(res_tmp.begin(), res_tmp.begin() + k);
		//out << "Success! Result num is: " << result.size() << endl;
	}
	else{
		//out << "Result num is: " << res_tmp.size() << endl;
		result = res_tmp;
		cout << "Not enough search results!" << endl;
	}

	query_result.resize(k);
	for (int i = 0; i < result.size(); i++){
		string img_path;
		if (i < 10)
			img_path = db_paths_cnn[result[i].second];
		
		else
			img_path = db_paths_fv[result[i].second];
			
		query_result[i] = make_pair(img_path, result[i].first);
	}

	clock_t end = clock();
	time += (end - start);

	// *************** Save Result ***************

	/*if (has_gt)
		showImages(query_point_cnn.first, result_path, result, groundtruth[splitFileName(query_point_cnn.first)], db_paths_cnn);
	else*/
		showImages(query_point_cnn.first, result_path, result, db_paths_cnn, db_paths_fv);


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