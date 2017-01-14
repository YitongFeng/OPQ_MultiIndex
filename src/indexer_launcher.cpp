// Copyright 2012 Yandex Artem Babenko
#include <iostream>
#include <fstream>

//#include <sys/time.h>
#include <time.h>
//#include <boost/program_options.hpp>

#include "indexer.h"
#include "perfomance_util.h"

//using namespace boost::program_options;
using std::ofstream;

int THREADS_COUNT;	// Number of threads for indexing

PointType point_type;	// Type, should be BVEC or FVEC

Dimensions SPACE_DIMENSION;	 // Number of coordinates in a point

string coarse_vocabs_file;  // File with vocabularies for multiindex structure

string fine_vocabs_file;	// File with vocabularies for reranking

string points_file;    // File with points to index
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
bool build_coarse_quantizations;
/**
* File with points coarse quantizations
*/
string coarse_quantizations_file;
/**
* How many points should we index
*/
int points_count;
/**
* Multiplicity of multiindex
*/
int multiplicity;

string result_path;

string index_param;

int SetOptions() {
	string data_name = "CNN2w";
	SPACE_DIMENSION = 256;
	points_count = 225533;

	string root_path = "I:\\Hashing\\Code\\output\\";	
	THREADS_COUNT = 1;
	multiplicity = 2;
	string index_path = root_path + data_name + "\\";	
	string result_path = root_path + data_name + "\\results\\";	
	files_prefix = index_path;
	points_file = index_path + data_name + "_base_NP.fvecs";	   // "index/fvImgDb/fgImgDb_base_NP.fvecs"
	metainfo_file = result_path + "metainfo.txt";
	coarse_vocabs_file = index_path + data_name + "_coarse.dat";    // "index/fvImgDb/fgImgDb_coarse.dat"
	fine_vocabs_file = index_path + data_name + "_fine.dat";
	coarse_quantizations_file = index_path + data_name + "_coarse_idx.dat";
	build_coarse_quantizations = 0;
	mode = USE_RESIDUALS;
	point_type = FVEC;

	return 0;
}

//float diff_timeval(timeval t1, timeval t2) {
//  return (float) (t1.tv_sec - t2.tv_sec) + (t1.tv_usec - t2.tv_usec) * 1e-6;
//}


int main() {

	SetOptions();

	cout << "Options are set ...\n";
	vector<Centroids> coarse_vocabs;
	vector<Centroids> fine_vocabs;
	ReadVocabularies<float>(coarse_vocabs_file, SPACE_DIMENSION, &coarse_vocabs);
	cout << "read coarse file\n";
	ReadFineVocabs<float>(fine_vocabs_file, &fine_vocabs);
	cout << "Vocs are read ...\n";
	
	clock_t start = clock();
	if (fine_vocabs.size() == 8) {
		MultiIndexer<RerankADC8> indexer(multiplicity);
		indexer.BuildMultiIndex(points_file, metainfo_file, points_count, coarse_vocabs,
			fine_vocabs, mode, build_coarse_quantizations,
			files_prefix, coarse_quantizations_file);
	}
	else if (fine_vocabs.size() == 16) {
		MultiIndexer<RerankADC16> indexer(multiplicity);
		indexer.BuildMultiIndex(points_file, metainfo_file, points_count, coarse_vocabs,
			fine_vocabs, mode, build_coarse_quantizations,
			files_prefix, coarse_quantizations_file);
	}
	float time = 0;
	clock_t end = clock();
	time += (end - start) / CLOCKS_PER_SEC / 60.0;

	ofstream out;
	out.open(result_path + "index_report.txt", ios::app);
	cout << " index_time " << time << "minutes" << endl;
	out << " index_time " << time << "minutes" << endl;
	return 0;
}
