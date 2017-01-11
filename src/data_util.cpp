// Copyright 2012 Yandex Artem Babenko

#include "data_util.h"

void splitString(string inString, vector<string> &sp){
	inString += " ";
	int size = inString.size();
	int st = 0;
	int pos = inString.find(" ", 0);
	while (pos<size && pos != inString.npos){
		string tmp = inString.substr(st, pos - st);
		st = pos + 1;
		//最后一个是空格
		if (st<size){
			pos = inString.find(" ", st);
		}
		else{
			break;
		}
		//防止最后有很多空格
		if (tmp != " " && tmp != ":"){
			sp.push_back(tmp);
		}
	}
}

void loadGt(const string gt_file, string imgFold, unordered_map<string, vector<pair<string, int>>>& gtOrder){
	// 加载 gt
	ifstream gtFile(gt_file);
	vector<string> gtList;
	int gtNum = 0;
	string tmpTline;
	while (getline(gtFile, tmpTline)){

		//gtFile >> tmpTline;
		if (tmpTline == "\n" || tmpTline == ""){ continue; }
		tmpTline += " ";
		vector<string > tmp;
		splitString(tmpTline, tmp);
		if (tmp.empty())	continue;

		/*for (int k = 0; k < tmp.size(); k++){
		tmp[k] = imgFold + "\\" + tmp[k];
		}*/

		string gtName = tmp[0];
		gtList.push_back(gtName);
		vector<pair<string, int>> tmpQ;
		for (int i = 0; i < tmp.size() - 1; i++){
			pair<string, int> tmpP;
			// 此gtfile含有原图 需要去除原图
			if (tmp[i + 1] != tmp[0]){
				tmpP = make_pair(tmp[i + 1], 1);
				tmpQ.push_back(tmpP);
			}
		}

		gtOrder[gtName] = tmpQ;


		gtNum++;
		//cout << "gen gt num: " << gtNum << "  : " << gtName << endl;
		//cout << gtOrder[gtName][0].first << endl;
	}
	//fclose(gtFile);
}

void loadDbOrQuery(const string data_file, const string pathfile, vector<pair<string, Point>>& data, PointType query_point_type){
	if (data.size() == 0)
		cerr << "DB size is 0, didn't initialize!" << endl;

	Points data_feas;
	if (query_point_type == BVEC) {
		ReadPoints<unsigned char, Coord>(data_file, &data_feas, data.size());
	}
	else if (query_point_type == FVEC) {
		ReadPoints<float, Coord>(data_file, &data_feas, data.size());
	}

	vector<string> paths;
	paths.resize(data.size());
	ifstream path_file(pathfile);
	for (int i = 0; i < paths.size(); i++)
		path_file >> paths[i];

	data.resize(paths.size());
	for (int i = 0; i < data.size(); i++){
		data[i] = make_pair(paths[i], data_feas[i]);
	}
	cout << "Data are read ...\n";
}

Distance Eucldistance(const Point& x, const Point& y) {
	Distance result = 0;
	Distance current_coord_diff;
	for (Dimensions d = 0; d < x.size(); ++d){
		current_coord_diff = x[d] - y[d];
		result += current_coord_diff * current_coord_diff;
	}
	return result;
}

Distance Eucldistance(const Point& x, const Point& y,
	const Dimensions start, const Dimensions finish) {
	Distance result = 0;
	Distance current_coord_diff;
	for (Dimensions d = start; d < finish; ++d){
		current_coord_diff = x[d] - y[d - start];
		result += current_coord_diff * current_coord_diff;
	}
	return result;
}

void GetSubpoints(const Points& points,
	const Dimensions start_dim,
	const Dimensions final_dim,
	Points* subpoints) {
	if (final_dim < start_dim) {
		throw std::logic_error("Final dim < Start dim");
	}
	subpoints->resize(points.size());
	for (PointId pid = 0; pid < points.size(); ++pid) {
		subpoints->at(pid).resize(final_dim - start_dim);
		for (Dimensions dim = start_dim; dim < final_dim; ++dim) {
			subpoints->at(pid)[dim] = points[pid][start_dim + dim];
		}
	}
}

ClusterId GetNearestClusterId(const Point& point,
	const Centroids& centroids,
	const Dimensions start_dim,
	const Dimensions final_dim) {
	if (final_dim < start_dim) {
		throw std::logic_error("Final dim < Start dim");
	}
	ClusterId nearest = 0;
	Distance min_distance = Eucldistance(point, centroids[0], start_dim, final_dim);
	for (PointId pid = 1; pid < centroids.size(); ++pid) {
		Distance current_distance = 0;
		current_distance = Eucldistance(point, centroids[pid], start_dim, final_dim);
		if (current_distance < min_distance) {
			min_distance = current_distance;
			nearest = pid;
		}
	}
	return nearest;
}

void GetResidual(const Point& point, const CoarseQuantization& coarse_quantizations,
	const vector<Centroids>& centroids, Point* residual)
{
	residual->resize(point.size());
	Dimensions subvector_dimension = point.size() / centroids.size();
	cblas_saxpy(point.size(), 1, &(point[0]), 1, &(residual->at(0)), 1);
	for (int subvector_index = 0; subvector_index < centroids.size(); ++subvector_index) {
		Dimensions start_dim = subvector_index * subvector_dimension;
		const Point& current_coarse_centroid = centroids[subvector_index][coarse_quantizations[subvector_index]];
		cblas_saxpy(subvector_dimension, -1, &(current_coarse_centroid[0]), 1, &(residual->at(start_dim)), 1);
	}
}

void GetResidual(const Point& point, const CoarseQuantization& coarse_quantizations,
	const vector<Centroids>& centroids, Coord* residual) {
	Dimensions subvector_dimension = point.size() / centroids.size();
	cblas_scopy(point.size(), &(point[0]), 1, residual, 1);
	for (int subvector_index = 0; subvector_index < centroids.size(); ++subvector_index) {
		Dimensions start_dim = subvector_index * subvector_dimension;
		const Point& current_coarse_centroid = centroids[subvector_index][coarse_quantizations[subvector_index]];
		cblas_saxpy(subvector_dimension, -1, &(current_coarse_centroid[0]), 1, &(residual[start_dim]), 1);
	}
}

void GetNearestClusterIdsForPointSubset(const Points& points, const Centroids& centroids,
	const PointId start_pid, const PointId final_pid,
	vector<ClusterId>* nearest) {
	if (final_pid < start_pid) {
		throw std::logic_error("Final pid < Start pid");
	}
	cout << start_pid << " point processing started\n";
	for (PointId pid = start_pid; pid < final_pid; ++pid) {
		if (pid % 10000 == 0) {
			cout << pid << endl;
		}
		nearest->at(pid) = GetNearestClusterId(points[pid], centroids, 0, points[0].size());
	}
	cout << final_pid << " point processing finished\n";
}

void GetNearestClusterIdsForSubpoints(const Points& points, const Centroids& centroids,
	const Dimensions start_dim, const Dimensions final_dim,
	int threads_count, vector<ClusterId>* nearest) {
	if (final_dim < start_dim) {
		throw std::logic_error("Final dim < Start dim");
	}
	cout << "Start getting nearest Cluster Ids..." << endl;
	Points subpoints;
	GetSubpoints(points, start_dim, final_dim, &subpoints);
	boost::thread_group threads;
	int subpoints_count = points.size() / threads_count;
	for (int thread_id = 0; thread_id < threads_count; ++thread_id) {
		PointId start_pid = subpoints_count * thread_id;
		PointId final_pid = start_pid + subpoints_count;
		threads.create_thread(boost::bind(&GetNearestClusterIdsForPointSubset, subpoints, centroids,
			start_pid, final_pid, nearest));
	}
	threads.join_all();
	cout << "Finish getting nearest Cluster Ids..." << endl;
}

void GetPointsCoarseQuaintizations(const Points& points,
	const vector<Centroids>& centroids,
	const int threads_count,
	vector<CoarseQuantization>* coarse_quantizations) {
	int number_of_subvectors = centroids.size();
	coarse_quantizations->resize(points.size(), CoarseQuantization(number_of_subvectors));
	Dimensions subvector_dimension = points[0].size() / number_of_subvectors;
	for (int centroids_index = 0; centroids_index < number_of_subvectors; ++centroids_index) {
		vector<ClusterId> cluster_labels;
		cluster_labels.resize(points.size());
		Dimensions start_dim = centroids_index * subvector_dimension;
		Dimensions final_dim = std::min((Dimensions)points[0].size(), start_dim + subvector_dimension);
		GetNearestClusterIdsForSubpoints(points, centroids[centroids_index],
			start_dim, final_dim, threads_count, &cluster_labels);
		for (PointId pid = 0; pid < points.size(); ++pid) {
			coarse_quantizations->at(pid)[centroids_index] = cluster_labels[pid];
		}
	}
}