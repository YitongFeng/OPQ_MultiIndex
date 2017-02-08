/** @file */
// Copyright 2012 Yandex Artem Babenko

//#pragma once
#ifndef SEARCHER_H_
#define SEARCHER_H_

#include <algorithm>
#include <map>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/vector.hpp>

//#include <mkl_cblas.h>
//#include "cblas.h"

#include "data_util.h"
#include "ordered_lists_merger.h"
//#include "perfomance_util.h"

//extern int THREADS_COUNT;

Dimensions SPACE_DIMENSION;

extern enum PointType point_type;

//vector<vector <float> > dq_distance;
//int query_id;

/**
* \typedef This typedef is used in the first stage of search when
* we get nearest centroids for each coarse subpace
*/
typedef vector<pair<Distance, ClusterId> > NearestSubspaceCentroids;

/**
* This is the main class for nearest neighbour search using multiindex
*/
template<class Record, class MetaInfo>
class MultiSearcher {
public:

	/**
	* Default constructor
	*/
	MultiSearcher();
	/**
	* Initiation function
	* @param index_files_prefix prefix of multiindex files providing the search
	* @param coarse_vocabs_filename file with coarse vocabs
	* @param fine_vocabs_filename file with fine vocabs for reranking
	* @param mode reranking approach
	* @param do_rerank should algorithm rerank short list or not
	*/
	/*void Init(const string& index_files_prefix,
		const string& coarse_vocabs_filename,
		const string& fine_vocabs_filename,
		const RerankMode& mode,
		const int subspace_centroids_to_consider,
		bool do_rerank, int use_originaldata);*/
	void Init(const string& index_files_prefix,
		const string& coarse_vocabs_filename_fv, const string& coarse_vocabs_filename_cnn,
		const string& fine_vocabs_filename_fv, const string& fine_vocabs_filename_cnn,
		const RerankMode& mode, const int subspace_centroids_to_consider,
		const bool do_rerank, int use_originaldata);
	
	void setmulti(const Multitable<int>& cell_edges_in, const vector<Record>& multiindex_in, 
		const vector<Centroids>& coarse_vocabs_in_, const vector<Centroids>& fine_vocabs_in_, 
		const vector<float*>& coarse_vocabs_matrices_in_, const vector<vector<float> >& coarse_centroids_norms_in_, int space_dim_in);
	/**
	* Main interface function
	* @param point query point
	* @param k number of neighbours to get
	* @param subpace_centroids_to_consider it defines the size of working index table
	* @param neighbours result - vector of point identifiers ordered by increasing of distance to query
	*/
	//void GetNearestNeighbours(const Point& point, int k, vector<pair<Distance, MetaInfo> >* neighbours) const;

	//add for real reranking
	/*void GetNearestNeighbours(const Point& point, int k,
		vector<pair<Distance, MetaInfo> >* neighbours, const Points& dataset) const;*/
	void GetNearestNeighbours(const Point& point, int k,
		vector<pair<Distance, MetaInfo> >* neighbours, const Points& dataset, int fea_flag);
	/**
	* Returns searcher perfomance tester
	*/
	//PerfTester& GetPerfTester();
	void setDebugTrue();
private:
	/**
	* This functions deserializes all structures for search
	* @param index_files_prefix prefix of multiindex files providing the search
	* @param coarse_vocabs_filename file with coarse vocabs
	* @param fine_vocabs_filename file with fine vocabs for reranking
	*/
	/*void DeserializeData(const string& index_files_prefix,
		const string& coarse_vocabs_filename,
		const string& fine_vocabs_filename);*/
	void DeserializeData(const string& index_files_prefix,
		const string& coarse_vocabs_filename_fv, const string& coarse_vocabs_filename_cnn,
		const string& fine_vocabs_filename_fv, const string& fine_vocabs_filename_cnn);
	/**
	* Function gets some nearest centroids for each coarse subspace
	* @param point query point
	* @param subspace_centroins_count how many nearest subcentroids to get
	* @param subspaces_short_lists result
	*/
	void GetNearestSubspacesCentroids(const Point& point,
		const int subspace_centroins_count,
		vector<NearestSubspaceCentroids>* subspaces_short_lists) const;

	/**
	* This fuctions traverses another cell of multiindex table
	* @param point query point
	* @param nearest_subpoints vector algorithm adds nearest neighbours in
	*/
	/*bool TraverseNextMultiIndexCell(const Point& point,
		vector<pair<Distance, MetaInfo> >* nearest_subpoints, const Points& dataset) const;*/
		bool TraverseNextMultiIndexCell(const Point& point, vector<pair<Distance, MetaInfo> >* nearest_subpoints, 
			const Points& dataset, OrderedListsMerger<Distance, ClusterId>& merger_, int& found_neghbours_count_) const;
	/**
	* This fuctions converts cells coordinates to appropriate range in array
	* @param cell_coordinates coordinates of the cell
	* @param cell_start first index of range
	* @param cell_finish last index of range
	*/
	inline void GetCellEdgesInMultiIndexArray(const vector<int>& cell_coordinates,
		int* cell_start, int* cell_finish) const;
	/**
	* This fuctions converts complex objects to arrays and
	* pointers for usage in BLAS
	*/
	void InitBlasStructures();
	/**
	* Lists of coarse centroids
	*/
	vector<Centroids> coarse_vocabs_;
	vector<Centroids> coarse_vocabs_fv_;
	vector<Centroids> coarse_vocabs_cnn_;
	/**
	* Lists of fine centroids
	*/
	vector<Centroids> fine_vocabs_;
	vector<Centroids> fine_vocabs_fv_;
	vector<Centroids> fine_vocabs_cnn_;
	/**
	* Merger for ordered merging subspaces centroids lists
	*/
	//mutable OrderedListsMerger<Distance, ClusterId> merger_;
	/**
	* Should algorithm use reranking or not
	*/
	bool do_rerank_;
	/**
	* Searcher perfomance tester
	*/
	//mutable PerfTester perf_tester_;
	/**
	* Common prefix of every index files
	*/
	string index_files_prefix_;
	/**
	* Multiindex data structures
	*/
	MultiIndex<Record> multiindex_;
	/**
	* Reranking approach
	*/
	RerankMode rerank_mode_;
	/**
	* Struct for BLAS
	*/
	vector<float*> coarse_vocabs_matrices_;
	vector<float*> coarse_vocabs_matrices_fv_;
	vector<float*> coarse_vocabs_matrices_cnn_;

	vector<vector<float> > coarse_centroids_norms_;
	vector<vector<float> > coarse_centroids_norms_fv_;
	vector<vector<float> > coarse_centroids_norms_cnn_;
	//mutable Coord* products_;
	//mutable vector<Coord> query_norms_;
	//mutable float* residual_;
	/**
	* Number of nearest to query centroids
	* to consider for each dimension
	*/
	int subspace_centroids_to_consider_;

	int use_originaldata_;

	int debug;	// 是否输出调试信息
	ofstream of_debug;
	//
	//mutable int query_id;
	/**
	* Number of neighbours found to this moment
	*/
	//mutable int found_neghbours_count_;
};
template<class Record, class MetaInfo>
void MultiSearcher<Record, MetaInfo>::setDebugTrue(){
	debug = 1;
	of_debug.open("search_debug.txt");
}

template<class Record, class MetaInfo>
inline void RecordToMetainfoAndDistance(const Coord* point,
	const Record& record,
	pair<Distance, MetaInfo>* result,
	const vector<int>& cell_coordinates,
	const vector<Centroids>& fine_vocabs, const Points& dataset, int use_originaldata_) {
}

/////////////// IMPLEMENTATION /////////////////////
template<class Record, class MetaInfo>
void MultiSearcher<Record, MetaInfo>::setmulti(const Multitable<int>& cell_edges_in, 
	const vector<Record>& multiindex_in, const vector<Centroids>& coarse_vocabs_in_, 
	const vector<Centroids>& fine_vocabs_in_, const vector<float*>& coarse_vocabs_matrices_in_, 
	const vector<vector<float> >& coarse_centroids_norms_in_, int space_dim_in){

	multiindex_.cell_edges = Multitable<int>(cell_edges_in.dimensions);
	multiindex_.multiindex.resize(multiindex_in.size());
	coarse_vocabs_.resize(coarse_vocabs_in_.size());
	fine_vocabs_.resize(fine_vocabs_in_.size());
	coarse_vocabs_matrices_.resize(coarse_vocabs_matrices_in_.size());
	coarse_centroids_norms_.resize(coarse_centroids_norms_in_.size());

	copy(cell_edges_in.table.begin(), cell_edges_in.table.end(), multiindex_.cell_edges.table.begin());
	copy(multiindex_in.begin(), multiindex_in.end(), multiindex_.multiindex.begin());
	copy(coarse_vocabs_in_.begin(), coarse_vocabs_in_.end(), coarse_vocabs_.begin());
	copy(fine_vocabs_in_.begin(), fine_vocabs_in_.end(), fine_vocabs_.begin());
	copy(coarse_vocabs_matrices_in_.begin(), coarse_vocabs_matrices_in_.end(), coarse_vocabs_matrices_.begin());
	copy(coarse_centroids_norms_in_.begin(), coarse_centroids_norms_in_.end(), coarse_centroids_norms_.begin());
	SPACE_DIMENSION = space_dim_in;
}

template<class Record, class MetaInfo>
MultiSearcher<Record, MetaInfo>::MultiSearcher() {
}

//template<class Record, class MetaInfo>
//void MultiSearcher<Record, MetaInfo>::DeserializeData(const string& index_files_prefix,
//	const string& coarse_vocabs_filename,
//	const string& fine_vocabs_filename) {
//	cout << "Data deserializing started...\n";
//
//	std::ifstream cell_edges(string(index_files_prefix + "_cell_edges.bin").c_str(), ios::binary);
//	if (!cell_edges.good() || !cell_edges || !cell_edges.is_open()) {
//		throw std::logic_error("Bad input cell edges stream");
//	}
//	boost::archive::binary_iarchive arc_cell_edges(cell_edges);
//	arc_cell_edges >> multiindex_.cell_edges;
//	cout << "Cell edges deserialized...\n";
//
//	ifstream multi_array(string(index_files_prefix + "_multi_array.bin").c_str(), ios::binary);
//	if (!multi_array.good() || !multi_array) {
//		throw std::logic_error("Bad input cell edges stream");
//	}
//	boost::archive::binary_iarchive arc_multi_array(multi_array);
//	arc_multi_array >> multiindex_.multiindex;
//	cout << "Multiindex deserialized...\n";
//	ReadVocabularies<float>(coarse_vocabs_filename, SPACE_DIMENSION, &coarse_vocabs_);
//	/*for(int i=0;i<2;i++)
//	{
//	int sum=0;
//	for(int j=0;j<coarse_vocabs_[i].size();j++)
//	{
//	sum+=coarse_vocabs_[i][j].size();
//	cout<<i<<" "<<j<<" "<<coarse_vocabs_[i][j].size()<<endl;
//	}
//	cout<<"sum:"<<sum<<endl;
//	}*/
//	cout << "Coarse vocabs deserialized...\n";
//	ReadFineVocabs<float>(fine_vocabs_filename, &fine_vocabs_);
//	cout << "Fine vocabs deserialized...\n";
//}

template<class Record, class MetaInfo>
void MultiSearcher<Record, MetaInfo>::DeserializeData(const string& index_files_prefix,
	const string& coarse_vocabs_filename_fv, const string& coarse_vocabs_filename_cnn,
	const string& fine_vocabs_filename_fv, const string& fine_vocabs_filename_cnn) {
	cout << "Data deserializing started...\n";

	std::ifstream if_cell_edges_fv(string(index_files_prefix + "_cell_edges_fv.bin").c_str(), ios::binary);
	if (!if_cell_edges_fv.good() || !if_cell_edges_fv || !if_cell_edges_fv.is_open()) {
		cout << "Bad input fv cell edges stream" << endl;
		throw std::logic_error("Bad input fv cell edges stream");
	}
	boost::archive::binary_iarchive arc_cell_edges_fv(if_cell_edges_fv);
	arc_cell_edges_fv >> multiindex_.cell_edges_fv;
	cout << "FV Cell edges deserialized...\n";

	std::ifstream if_cell_edges_cnn(string(index_files_prefix + "_cell_edges_cnn.bin").c_str(), ios::binary);
	if (!if_cell_edges_cnn.good() || !if_cell_edges_cnn || !if_cell_edges_cnn.is_open()) {
		cout << "Bad input cnn cell edges stream" << endl;
		throw std::logic_error("Bad input cnn cell edges stream");
	}
	boost::archive::binary_iarchive arc_cell_edges_cnn(if_cell_edges_cnn);
	arc_cell_edges_cnn >> multiindex_.cell_edges_cnn;
	cout << "CNN Cell edges deserialized...\n";

	ifstream multi_array_fv(string(index_files_prefix + "_multi_array_fv.bin").c_str(), ios::binary);
	if (!multi_array_fv.good() || !multi_array_fv) {
		throw std::logic_error("Bad input multi_array stream");
	}
	boost::archive::binary_iarchive arc_multi_array_fv(multi_array_fv);
	arc_multi_array_fv >> multiindex_.multiindex_fv;
	cout << "FV Multiindex deserialized...\n";
	ReadVocabularies<float>(coarse_vocabs_filename_fv, FV_DIMENSION, &coarse_vocabs_fv_);

	ifstream multi_array_cnn(string(index_files_prefix + "_multi_array_cnn.bin").c_str(), ios::binary);
	if (!multi_array_cnn.good() || !multi_array_cnn) {
		throw std::logic_error("Bad input cell edges stream");
	}
	boost::archive::binary_iarchive arc_multi_array_cnn(multi_array_cnn);
	arc_multi_array_cnn >> multiindex_.multiindex_cnn;
	cout << "CNN Multiindex deserialized...\n";
	ReadVocabularies<float>(coarse_vocabs_filename_cnn, CNN_DIMENSION, &coarse_vocabs_cnn_);
	/*for(int i=0;i<2;i++)
	{
	int sum=0;
	for(int j=0;j<coarse_vocabs_[i].size();j++)
	{
	sum+=coarse_vocabs_[i][j].size();
	cout<<i<<" "<<j<<" "<<coarse_vocabs_[i][j].size()<<endl;
	}
	cout<<"sum:"<<sum<<endl;
	}*/
	cout << "Coarse vocabs deserialized...\n";
	ReadFineVocabs<float>(fine_vocabs_filename_fv, &fine_vocabs_fv_);
	ReadFineVocabs<float>(fine_vocabs_filename_cnn, &fine_vocabs_cnn_);
	cout << "Fine vocabs deserialized...\n";
}

//template<class Record, class MetaInfo>
//void MultiSearcher<Record, MetaInfo>::Init(const string& index_files_prefix,
//	const string& coarse_vocabs_filename,
//	const string& fine_vocabs_filename,
//	const RerankMode& mode,
//	const int subspace_centroids_to_consider,
//	const bool do_rerank, int use_originaldata) {
//	do_rerank_ = do_rerank;
//	use_originaldata_ = use_originaldata;
//	index_files_prefix_ = index_files_prefix;
//	subspace_centroids_to_consider_ = subspace_centroids_to_consider;
//	DeserializeData(index_files_prefix, coarse_vocabs_filename, fine_vocabs_filename);
//	rerank_mode_ = mode;
//	//merger_.GetYieldedItems().table.resize(std::pow((float)subspace_centroids_to_consider,
//	//	(int)coarse_vocabs_.size()));	// multi table resize to (2*k)^2 eg. for 10-NN is 20^2=400
//	//for (int i = 0; i < coarse_vocabs_.size(); ++i) {
//	//	merger_.GetYieldedItems().dimensions.push_back(subspace_centroids_to_consider);
//	//}
//	InitBlasStructures();
//}

template<class Record, class MetaInfo>
void MultiSearcher<Record, MetaInfo>::Init(const string& index_files_prefix,
	const string& coarse_vocabs_filename_fv, const string& coarse_vocabs_filename_cnn,
	const string& fine_vocabs_filename_fv, const string& fine_vocabs_filename_cnn,
	const RerankMode& mode, const int subspace_centroids_to_consider,
	const bool do_rerank, int use_originaldata) {
	debug = 0;
	do_rerank_ = do_rerank;
	use_originaldata_ = use_originaldata;
	index_files_prefix_ = index_files_prefix;
	subspace_centroids_to_consider_ = subspace_centroids_to_consider;
	DeserializeData(index_files_prefix, coarse_vocabs_filename_fv, coarse_vocabs_filename_cnn, 
		fine_vocabs_filename_fv, fine_vocabs_filename_cnn);
	rerank_mode_ = mode;
	//merger_.GetYieldedItems().table.resize(std::pow((float)subspace_centroids_to_consider,
	//	(int)coarse_vocabs_.size()));	// multi table resize to (2*k)^2 eg. for 10-NN is 20^2=400
	//for (int i = 0; i < coarse_vocabs_.size(); ++i) {
	//	merger_.GetYieldedItems().dimensions.push_back(subspace_centroids_to_consider);
	//}
	InitBlasStructures();
}

template<class Record, class MetaInfo>
void MultiSearcher<Record, MetaInfo>::InitBlasStructures(){
	coarse_vocabs_matrices_fv_.resize(coarse_vocabs_fv_.size());
	coarse_vocabs_matrices_cnn_.resize(coarse_vocabs_cnn_.size());
	coarse_centroids_norms_fv_.resize(coarse_vocabs_fv_.size(), vector<float>(coarse_vocabs_fv_[0].size()));
	coarse_centroids_norms_cnn_.resize(coarse_vocabs_cnn_.size(), vector<float>(coarse_vocabs_cnn_[0].size()));

	for (int coarse_id = 0; coarse_id < coarse_vocabs_matrices_fv_.size(); ++coarse_id) {
		coarse_vocabs_matrices_fv_[coarse_id] = new float[coarse_vocabs_fv_[0].size() * coarse_vocabs_fv_[0][0].size()];
		for (int i = 0; i < coarse_vocabs_fv_[0].size(); ++i) {
			Coord norm = 0;
			for (int j = 0; j < coarse_vocabs_fv_[0][0].size(); ++j) {
				coarse_vocabs_matrices_fv_[coarse_id][coarse_vocabs_fv_[0][0].size() * i + j] = coarse_vocabs_fv_[coarse_id][i][j];
				norm += coarse_vocabs_fv_[coarse_id][i][j] * coarse_vocabs_fv_[coarse_id][i][j];
			}
			coarse_centroids_norms_fv_[coarse_id][i] = norm;
		}
	}

	for (int coarse_id = 0; coarse_id < coarse_vocabs_matrices_cnn_.size(); ++coarse_id) {
		coarse_vocabs_matrices_cnn_[coarse_id] = new float[coarse_vocabs_cnn_[0].size() * coarse_vocabs_cnn_[0][0].size()];
		for (int i = 0; i < coarse_vocabs_cnn_[0].size(); ++i) {
			Coord norm = 0;
			for (int j = 0; j < coarse_vocabs_cnn_[0][0].size(); ++j) {
				coarse_vocabs_matrices_cnn_[coarse_id][coarse_vocabs_cnn_[0][0].size() * i + j] = coarse_vocabs_cnn_[coarse_id][i][j];
				norm += coarse_vocabs_cnn_[coarse_id][i][j] * coarse_vocabs_cnn_[coarse_id][i][j];
			}
			coarse_centroids_norms_cnn_[coarse_id][i] = norm;
		}
	}
	//Coord* products_ = new Coord[coarse_vocabs_[0].size()];
	//vector<Coord> query_norms_.resize(coarse_vocabs_[0].size());
	//float* residual_ = new Coord[coarse_vocabs_[0][0].size() * coarse_vocabs_.size()];
}

//template<class Record, class MetaInfo>
//PerfTester& MultiSearcher<Record, MetaInfo>::GetPerfTester() {
//	return perf_tester_;
//}

template<class Record, class MetaInfo>
void MultiSearcher<Record, MetaInfo>::GetNearestSubspacesCentroids(const Point& point,
	const int subspace_centroins_count,
	vector<NearestSubspaceCentroids>*
	subspaces_short_lists) const {
	std::stringstream aa;
	subspaces_short_lists->resize(coarse_vocabs_.size());
	Dimensions subspace_dimension = point.size() / coarse_vocabs_.size();
	for (int subspace_index = 0; subspace_index < coarse_vocabs_.size(); ++subspace_index) {
		Dimensions start_dim = subspace_index * subspace_dimension;
		Dimensions final_dim = min((Dimensions)point.size(), start_dim + subspace_dimension);
		Coord query_norm = cblas_sdot(final_dim - start_dim, &(point[start_dim]), 1, &(point[start_dim]), 1);

		vector<Coord> query_norms_(coarse_vocabs_[0].size());
		std::fill(query_norms_.begin(), query_norms_.end(), query_norm);
		cblas_saxpy(coarse_vocabs_[0].size(), 1, &(coarse_centroids_norms_[subspace_index][0]), 1, &(query_norms_[0]), 1);
		cout << "coarse vocab matrices size: " << coarse_vocabs_matrices_.size() << endl;
		cout << "sbspace index : " << subspace_index << endl;
		cblas_sgemv(CblasRowMajor, CblasNoTrans, coarse_vocabs_[0].size(), subspace_dimension, -2.0,
			coarse_vocabs_matrices_[subspace_index], subspace_dimension, &(point[start_dim]), 1, 1, &(query_norms_[0]), 1);
		subspaces_short_lists->at(subspace_index).resize(query_norms_.size());
		for (int i = 0; i < query_norms_.size(); ++i) {
			subspaces_short_lists->at(subspace_index)[i] = std::make_pair(query_norms_[i], i);
		}
		std::nth_element(subspaces_short_lists->at(subspace_index).begin(),
			subspaces_short_lists->at(subspace_index).begin() + subspace_centroins_count,
			subspaces_short_lists->at(subspace_index).end());
		subspaces_short_lists->at(subspace_index).resize(subspace_centroins_count);
		std::sort(subspaces_short_lists->at(subspace_index).begin(),
			subspaces_short_lists->at(subspace_index).end());
	}
}

template<class Record, class MetaInfo>
void MultiSearcher<Record, MetaInfo>::GetCellEdgesInMultiIndexArray(const vector<int>& cell_coordinates,
	int* cell_start, int* cell_finish) const {
	int global_index = multiindex_.cell_edges.GetCellGlobalIndex(cell_coordinates);
	*cell_start = multiindex_.cell_edges.table[global_index];
	if (global_index + 1 == multiindex_.cell_edges.table.size()) {
		*cell_finish = multiindex_.multiindex.size();
	}
	else {
		*cell_finish = multiindex_.cell_edges.table[global_index + 1];
	}
}


template<class Record, class MetaInfo>
bool MultiSearcher<Record, MetaInfo>::TraverseNextMultiIndexCell(const Point& point,vector<pair<Distance, MetaInfo> >* nearest_subpoints,
	const Points& dataset, OrderedListsMerger<Distance, ClusterId>& merger_, int& found_neghbours_count_) const {
	MergedItemIndices cell_inner_indices;
	clock_t before = clock();

	if (!merger_.GetNextMergedItemIndices(&cell_inner_indices)) {
		return false;
	}
	clock_t after = clock();
	//perf_tester_.cell_coordinates_time += after - before;
	vector<int> cell_coordinates(cell_inner_indices.size());	// save coarse quantization. int[2]
	for (int list_index = 0; list_index < merger_.lists_ptr->size(); ++list_index) {
		cell_coordinates[list_index] = merger_.lists_ptr->at(list_index)[cell_inner_indices[list_index]].second;
	}
	int cell_start, cell_finish;
	before = clock();
	GetCellEdgesInMultiIndexArray(cell_coordinates, &cell_start, &cell_finish);
	after = clock();
	//perf_tester_.cell_edges_time += after - before;
	if (cell_start >= cell_finish) {
		return true;
	}
	typename vector<Record>::const_iterator it = multiindex_.multiindex.begin() + cell_start;

	float* residual_ = new Coord[coarse_vocabs_[0][0].size() * coarse_vocabs_.size()];
	GetResidual(point, cell_coordinates, coarse_vocabs_, residual_);
	cell_finish = min((int)cell_finish, cell_start + (int)nearest_subpoints->size() - found_neghbours_count_);
	for (int array_index = cell_start; array_index < cell_finish; ++array_index) {
		if (rerank_mode_ == USE_RESIDUALS) {
			RecordToMetainfoAndDistance<Record, MetaInfo>(residual_, *it,
				&(nearest_subpoints->at(found_neghbours_count_)),
				cell_coordinates, fine_vocabs_, dataset, use_originaldata_);
		}
		else if (rerank_mode_ == USE_INIT_POINTS) {
			RecordToMetainfoAndDistance<Record, MetaInfo>(&(point[0]), *it,
				&(nearest_subpoints->at(found_neghbours_count_)),
				cell_coordinates, fine_vocabs_, dataset, use_originaldata_);
		}
		//perf_tester_.NextNeighbour();
		++found_neghbours_count_;
		++it;
	}
	return true;
}

//template<class Record, class MetaInfo>
//void MultiSearcher<Record, MetaInfo>::GetNearestNeighbours(const Point& point, int k,
//	vector<pair<Distance, MetaInfo> >* neighbours, const Points& dataset) const {
//
//	//perf_tester_.handled_queries_count += 1;
//	//perf_tester_.ResetQuerywiseStatistic();
//	//clock_t start = clock();
//	//perf_tester_.search_start = start;
//	//clock_t before = clock();
//
//	assert(k > 0);
//	neighbours->resize(k);
//	vector<NearestSubspaceCentroids> subspaces_short_lists;
//	assert(subspace_centroids_to_consider_ > 0);
//	GetNearestSubspacesCentroids(point, subspace_centroids_to_consider_, &subspaces_short_lists);  // ********** key ***********
//	cout << "get nearest subspace success! " << endl;
//	//clock_t after = clock();
//	//perf_tester_.nearest_subcentroids_time += after - before;
//
//	//clock_t before_merger = clock();
//	OrderedListsMerger<Distance, ClusterId> merger_(subspace_centroids_to_consider_, coarse_vocabs_.size());
//	merger_.setLists(subspaces_short_lists);
//	cout << "merge success! " << endl;
//	//clock_t after_merger = clock();
//	//perf_tester_.merger_init_time += after_merger - before_merger;
//
//	clock_t before_traversal = clock();
//	int found_neghbours_count_ = 0;
//	bool traverse_next_cell = true;
//	int cells_visited = 0;
//	while (found_neghbours_count_ < k && traverse_next_cell) {
//		//perf_tester_.cells_traversed += 1;
//		traverse_next_cell = TraverseNextMultiIndexCell(point, neighbours, dataset, merger_, found_neghbours_count_);	    // ******** key **********
//		cells_visited += 1;
//	}
//	if (found_neghbours_count_ < k){
//		std::sort(neighbours->begin(), neighbours->end(), [](pair<Distance, MetaInfo>& n1, pair<Distance, MetaInfo>& n2){
//			return n1.first > n2.first;
//		});
//		neighbours->resize(found_neghbours_count_);
//	}
//
//
//	clock_t after_traversal = clock();
//	//perf_tester_.full_traversal_time += after_traversal - before_traversal;
//
//	if (do_rerank_) {
//		std::sort(neighbours->begin(), neighbours->end());
//	}
//	clock_t finish = clock();
//	//perf_tester_.full_search_time += finish - start;
//
//}

template<class Record, class MetaInfo>
void MultiSearcher<Record, MetaInfo>::GetNearestNeighbours(const Point& point, int k,
	vector<pair<Distance, MetaInfo> >* neighbours, const Points& dataset, int fea_flag){
	// fea_flag cnn is 0, fv is 1
	if (fea_flag == 0){
		setmulti(multiindex_.cell_edges_cnn, multiindex_.multiindex_cnn, coarse_vocabs_cnn_, fine_vocabs_cnn_,
			coarse_vocabs_matrices_cnn_, coarse_centroids_norms_cnn_, CNN_DIMENSION);
	}
	else {
		setmulti(multiindex_.cell_edges_fv, multiindex_.multiindex_fv, coarse_vocabs_fv_, fine_vocabs_fv_,
			coarse_vocabs_matrices_fv_, coarse_centroids_norms_fv_, FV_DIMENSION);
	}

	//perf_tester_.handled_queries_count += 1;
	//perf_tester_.ResetQuerywiseStatistic();
	//clock_t start = clock();
	//perf_tester_.search_start = start;
	//clock_t before = clock();

	assert(k > 0);
	neighbours->resize(k);
	vector<NearestSubspaceCentroids> subspaces_short_lists;
	assert(subspace_centroids_to_consider_ > 0);
	GetNearestSubspacesCentroids(point, subspace_centroids_to_consider_, &subspaces_short_lists);  // ********** key ***********
	cout << "get nearest subspace success! " << endl;

	if(debug){
		of_debug << "fea flag: " << fea_flag << endl;
		of_debug << "query point: " << endl;
		for (int ii = 0; ii < 10; ii++){
			of_debug << point[ii] << " ";
		}
		of_debug << "\n nearest subspace: " << endl;
		of_debug << subspaces_short_lists.size() << endl;
		for (int ii = 0; ii < subspaces_short_lists.size(); ii++){
			of_debug << subspaces_short_lists[ii][0].second << " ";
		}
		of_debug << endl;
	}
	//clock_t after = clock();
	//perf_tester_.nearest_subcentroids_time += after - before;

	//clock_t before_merger = clock();
	OrderedListsMerger<Distance, ClusterId> merger_(subspace_centroids_to_consider_, coarse_vocabs_.size());
	merger_.setLists(subspaces_short_lists);
	cout << "merge success! " << endl;
	//clock_t after_merger = clock();
	//perf_tester_.merger_init_time += after_merger - before_merger;

	clock_t before_traversal = clock();
	int found_neghbours_count_ = 0;
	bool traverse_next_cell = true;
	int cells_visited = 0;
	while (found_neghbours_count_ < k && traverse_next_cell) {
		//perf_tester_.cells_traversed += 1;
		traverse_next_cell = TraverseNextMultiIndexCell(point, neighbours, dataset, merger_, found_neghbours_count_);	    // ******** key **********
		cells_visited += 1;
	}
	if (found_neghbours_count_ < k){
		std::sort(neighbours->begin(), neighbours->end(), [](pair<Distance, MetaInfo>& n1, pair<Distance, MetaInfo>& n2){
			return n1.first > n2.first;
		});
		neighbours->resize(found_neghbours_count_);
	}


	clock_t after_traversal = clock();
	//perf_tester_.full_traversal_time += after_traversal - before_traversal;

	if (do_rerank_) {
		std::sort(neighbours->begin(), neighbours->end());
	}
	clock_t finish = clock();
	//perf_tester_.full_search_time += finish - start;

}

template<>
inline void RecordToMetainfoAndDistance<RerankADC8, PointId>(const Coord* point, const RerankADC8& record,
	pair<Distance, PointId>* result,
	const vector<int>& cell_coordinates,
	const vector<Centroids>& fine_vocabs, const Points& dataset, int use_originaldata_) {
	
	assert(SPACE_DIMENSION > 0);
	result->second = record.pid;
	if (use_originaldata_ == 1)
	{
		Distance distance = 0;
		Point data = dataset[record.pid];
		for (int i = 0; i<SPACE_DIMENSION; i++)
		{
			Coord diff = data[i] - point[i];
			distance += diff * diff;
		}
		result->first = distance;
	}
	else
	{
		int coarse_clusters_count = cell_coordinates.size();
		int fine_clusters_count = fine_vocabs.size();
		int coarse_to_fine_ratio = fine_clusters_count / coarse_clusters_count;
		int subvectors_dim = SPACE_DIMENSION / fine_clusters_count;
		char* rerank_info_ptr = (char*)&record + sizeof(record.pid);
		for (int centroid_index = 0; centroid_index < fine_clusters_count; ++centroid_index) {
			int start_dim = centroid_index * subvectors_dim;
			int final_dim = start_dim + subvectors_dim;
			FineClusterId pid_nearest_centroid = *((FineClusterId*)rerank_info_ptr);
			rerank_info_ptr += sizeof(FineClusterId);
			int current_coarse_index = centroid_index / coarse_to_fine_ratio;
			Distance subvector_distance = 0;
			for (int i = start_dim; i < final_dim; ++i) {
				Coord diff = fine_vocabs[centroid_index][pid_nearest_centroid][i - start_dim] - point[i];
				subvector_distance += diff * diff;
			}
			result->first += subvector_distance;
		}
	}

}

template<>
inline void RecordToMetainfoAndDistance<RerankADC16, PointId>(const Coord* point, const RerankADC16& record,
	pair<Distance, PointId>* result,
	const vector<int>& cell_coordinates,
	const vector<Centroids>& fine_vocabs, const Points& dataset, int use_originaldata_) {
	result->second = record.pid;
	if (use_originaldata_ == 1)
	{
		Distance distance = 0;
		Point data = dataset[record.pid];
		for (int i = 0; i<SPACE_DIMENSION; i++)
		{
			Coord diff = data[i] - point[i];
			distance += diff * diff;
		}
		result->first = distance;
	}
	else
	{
		int coarse_clusters_count = cell_coordinates.size();
		int fine_clusters_count = fine_vocabs.size();
		int coarse_to_fine_ratio = fine_clusters_count / coarse_clusters_count;
		int subvectors_dim = SPACE_DIMENSION / fine_clusters_count;
		char* rerank_info_ptr = (char*)&record + sizeof(record.pid);
		for (int centroid_index = 0; centroid_index < fine_clusters_count; ++centroid_index) {
			int start_dim = centroid_index * subvectors_dim;
			int final_dim = start_dim + subvectors_dim;
			FineClusterId pid_nearest_centroid = *((FineClusterId*)rerank_info_ptr);
			rerank_info_ptr += sizeof(FineClusterId);
			int current_coarse_index = centroid_index / coarse_to_fine_ratio;
			Distance subvector_distance = 0;
			for (int i = start_dim; i < final_dim; ++i) {
				Coord diff = fine_vocabs[centroid_index][pid_nearest_centroid][i - start_dim] - point[i];
				subvector_distance += diff * diff;
			}
			result->first += subvector_distance;
		}
	}

}

template class MultiSearcher<RerankADC8, PointId>;
template class MultiSearcher<RerankADC16, PointId>;
template class MultiSearcher<PointId, PointId>;

#endif

