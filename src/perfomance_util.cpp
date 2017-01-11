// Copyright 2012 Yandex Artem Babenko
#include "perfomance_util.h"

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

	for (int i = 0; i < 21; ++i) {
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
	if (current_points_count >= list_length_thresholds_[current_threshold_index_]) {
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
		<< (double)nearest_subcentroids_time / handled_queries_count / CLOCKS_PER_SEC << endl;
	out << "Average cache init time: "
		<< (double)cache_init_time / handled_queries_count / CLOCKS_PER_SEC << endl;
	out << "Average merger init time: "
		<< (double)merger_init_time / handled_queries_count / CLOCKS_PER_SEC << endl;
	out << "Average full traversal time: "
		<< (double)full_traversal_time / handled_queries_count / CLOCKS_PER_SEC << endl;
	out << "Average cells coordinates getting time: "
		<< (double)cell_coordinates_time / handled_queries_count / CLOCKS_PER_SEC << endl;
	out << "Average cell edges getting time: "
		<< (double)cell_edges_time / handled_queries_count / CLOCKS_PER_SEC << endl;
	out << "Average residual time: "
		<< (double)residual_time / handled_queries_count / CLOCKS_PER_SEC << endl;
	out << "Average refining time: "
		<< (double)refining_time / handled_queries_count / CLOCKS_PER_SEC << endl;
	out << "Average full search time: "
		<< (double)full_search_time / handled_queries_count / CLOCKS_PER_SEC << endl;
}

void PerfTester::DoReport() {
	std::ofstream out(report_file_.c_str(), ios::app);
	DoReport(out);
}

void removeDublicates(vector<DistanceToPoint>& res_tmp, vector<string>& db_paths){
	sort(res_tmp.begin(), res_tmp.end(), comp_with_pointid);
	for (int ii = 1; ii < res_tmp.size(); ii++){
		int id1 = res_tmp[ii - 1].second;
		int id2 = res_tmp[ii].second;
		if ((id2 == (id1 + 1) || id2 == id1) && db_paths[id1] == db_paths[id2]){	//只有pointid相邻的特征才可能是同一幅图的特征
			if (res_tmp[ii - 1].first < res_tmp[ii].first){
				res_tmp.erase(res_tmp.begin() + ii);
				ii--;
			}
			else{
				res_tmp.erase(res_tmp.begin() + ii - 1);
				ii--;
			}
		}
	}
}

void showImages(ImagePath queryName, string saveFold, vector< DistanceToPoint> qRe, vector<pair<string, int> > qOrder, vector<string> paths){
	if (qOrder.empty()){
		cerr << "ShowImage functon, gt is empty! " << queryName << endl;
	}

	ofstream outfile;
	outfile.open("result.txt", ios::app);
	string imgName = splitFileName(queryName);
	string imgSaveFold = saveFold + "\\" + imgName;
	string or = "md " + imgSaveFold;
	if (access(imgSaveFold.c_str(), 0) != 0){
		system(or.c_str());
	}

	// save query image
	string qSavePath = imgSaveFold + "\\" + "0_0_" + imgName;
	Mat imQ = cv::imread(queryName);
	/*bool Qbg = 0;
	vector<Rect> Qrect;
	proposal(&IplImage(imQ), Qrect, Qbg);
	for (int i = 0; i < Qrect.size(); i++){
	Rect rectShow = Qrect[i];
	rectangle(imQ, cvPoint(rectShow.x, rectShow.y), cvPoint(rectShow.x + rectShow.width, rectShow.y + rectShow.height), CV_RGB(255, 0, 0), 3, 8, 0);
	}*/
	imwrite(qSavePath, imQ);
	//copyImg(queryName, qSavePath);
	outfile << "query: " << imgName << endl;


	//保存label 图片
	for (int i = 0; i < qOrder.size(); i++){
		int labelOrder = 0;
		float disF = 0;
		int proposalOrder = 0;
		string a = qOrder[i].first;

		for (int j = 0; j < qRe.size(); j++){
			string b = splitFileName(paths[qRe[j].second]);
			if (a == b){
				labelOrder = j + 1;
				disF = qRe[j].first;
				proposalOrder = 0;
				break;
			}
		}
		char num[100], dis[100];
		sprintf(num, "%d", labelOrder);
		sprintf(dis, "%3f", disF);
		string labelNum(num);
		string disFS(dis);
		string name = a;
		string lSavePath = imgSaveFold + "\\" + "0_1_" + labelNum + "_" + disFS + "_" + name;


		cout << "label : " << name << " order: " << labelOrder << " dis : " << disF << endl;
		outfile << "label : " << name << " order: " << labelOrder << " dis : " << disF << endl;

		//为了显示的时候加上候选框		----2016-7-13
		/*Mat imTmp = imread(labelName);
		bool bg = 0;
		vector<Rect> rect;
		proposal(&IplImage(imTmp), rect, bg);
		Rect rectShow = rect[proposalOrder];
		rectangle(imTmp, cvPoint(rectShow.x, rectShow.y), cvPoint(rectShow.x + rectShow.width, rectShow.y + rectShow.height), CV_RGB(255, 0, 0), 3, 8, 0);

		imwrite(lSavePath, imTmp);*/
		//copyImg(labelName, lSavePath);
	}

	for (int i = 0; i < 20; i++){
		string path = paths[qRe[i].second];
		float disF = qRe[i].first;
		string name = splitFileName(path);
		int proposalOrder = 0;

		char num[100], dis[100];
		sprintf(num, "%d", i + 1);
		sprintf(dis, "%3f", disF);
		string numS(num), disS(dis);
		string rSavePath = imgSaveFold + "\\" + "0_2_" + numS + "_" + disS + "_" + name;

		cout << path << endl;
		cout << "result : " << i + 1 << " name :" << name << " dis : " << disF << endl;
		outfile << "result : " << i + 1 << " name :" << name << " dis : " << disF << endl;

		Mat imTmp = cv::imread(path);
		bool bg = 0;
		/*vector<Rect> rect;
		proposal(&IplImage(imTmp), rect, bg);
		Rect rectShow = rect[proposalOrder];
		rectangle(imTmp, cvPoint(rectShow.x, rectShow.y), cvPoint(rectShow.x + rectShow.width, rectShow.y + rectShow.height), CV_RGB(255, 0, 0), 3, 8, 0);
		*/
		imwrite(rSavePath, imTmp);

	}
	outfile.close();
}


int GetRecallAt(const int length, const vector<PointId>& groundtruth,
	const vector<DistanceToPoint>& result) {
	if (groundtruth.empty()) {
		cout << "Groundtruth is empty!" << endl;
		return 0;
	}
	for (int index = 0; index < length && index < result.size(); ++index) {
		if (result[index].second == groundtruth[0]) {
			return 1;
		}
	}
	return 0;
}

double GetPresicionAt(const int length, const vector<PointId>& groundtruth,
	const vector<DistanceToPoint>& result) {

	int found = 0;
	std::set<PointId> ground_points;
	for (int i = 0; i < groundtruth.size(); ++i) {
		ground_points.insert(groundtruth[i]);
	}
	for (int index = 0; index < length && index < result.size(); ++index) {
		if (ground_points.find(result[index].second) != ground_points.end()) {
			found += 1;
		}
	}
	return (double)found / length; //
}

double GetRecall(int k, const vector<PointId>& groundtruth,
	const vector<DistanceToPoint>& result) {
	if (groundtruth.empty()) {
		cout << "Groundtruth is empty!" << endl;
		return 0;
	}
	std::set<PointId> returned_points;
	for (int i = 0; i < k; ++i) {
		returned_points.insert(result[i].second);
	}
	double found = 0.0;
	for (int index = 0; index < k; ++index) {
		if (returned_points.find(groundtruth[index]) != returned_points.end()) {
			found += 1;
		}
	}
	return found / k;
}

double GetRecall(int k, vector<pair<string, int>>& gt_i,
	vector<DistanceToPoint>& res_i, vector<string>& db_paths) {
	if (gt_i.empty()) {
		cout << "Groundtruth is empty!" << endl;
		return 0;
	}
	std::set<string> returned_points_imgname;
	int i = 0;
	while (returned_points_imgname.size() < k){
		if (i >= res_i.size())
			cerr << "Index i of result out of range!" << endl;

		string imgname = splitFileName(db_paths[res_i[i].second]);
		returned_points_imgname.insert(imgname);
		i++;
	}
	double found = 0.0;
	assert(returned_points_imgname.size() == k);
	for (int index = 0; index < gt_i.size(); ++index) {
		if (returned_points_imgname.find(gt_i[index].first) != returned_points_imgname.end()) {
			found += 1;
		}
	}
	return found / gt_i.size();
}

float get_distance(Point data, Point query)
{
	float distance = 0.0;
	for (int i = 0; i< data.size(); i++)
	{
		float diff = data[i] - query[i];
		distance += diff * diff;
	}
	return distance;
}
float compute_relative_distance_error(int k, const Points& dataset, const Point& query, const vector<PointId>& groundtruth, const vector<DistanceToPoint>& result)
{
	float dist = 0.0;
	for (int i = 0; i <k; ++i) //dists.rows
	{
		float d_qr = get_distance(dataset[result[i].second], query);   //query.cols
		float d_qg = get_distance(dataset[groundtruth[i]], query);
		float d = (d_qr - d_qg) / d_qg;
		if (d <= 4)
			dist += d;
		else
			dist += 4;
	}
	return dist / k;
}

float compute_number_closer(int k, const vector<PointId>& groundtruth, const vector<DistanceToPoint>& result)
{
	float rate = 0.0;
	for (int gs_n = 0; gs_n < k; gs_n++)
	{
		for (int re_n = 0; re_n < k; re_n++)
		{
			if (groundtruth[gs_n] == result[re_n].second)
			{
				rate += (float)(gs_n + 1) / (re_n + 1);
				break;
			}
		}
	}
	return rate / k;
}

float compute_mean_reciprocal_rank(int k, const vector<PointId>& groundtruth, const vector<DistanceToPoint>& result)
{
	float sum = 0;
	float rate = 0.0;
	for (int re_n = 0; re_n < k; re_n++)
	{
		if (groundtruth[0] == result[re_n].second)
		{
			rate += 1.0 / (re_n + 1);
			break;
		}
	}
	return rate;
}

float compute_mean_average_precision(int k, const vector<PointId>& groundtruth, const vector<DistanceToPoint>& result)
{
	float rate = 0.0;
	int found_last = 0;
	std::set<PointId> gnds;
	for (int i = 0; i <k; i++)
	{
		gnds.insert(groundtruth[i]);
		int count = 0;
		for (int j = 0; j <= i; j++)
			if (gnds.find(result[j].second) != gnds.end())
			{
				count++;
			}
		rate += count*1.0 * (count - found_last) / (i + 1);
		found_last = count;
	}
	return rate / k;
}

float compute_discounted_culmulative_gain(int k, const vector<PointId>& groundtruth, const vector<DistanceToPoint>& result)
{
	float rate = 0.0;

	std::set<PointId> gnds;
	for (int i = 0; i <k; ++i) {
		gnds.insert(groundtruth[i]);
	}
	for (unsigned i = 0; i <k; i++)
	{
		if (gnds.find(result[i].second) != gnds.end())
		{
			rate += 1.0 / log2(i + 2);
		}
	}
	return rate;
}

