#include <fstream>
#include <vector>
// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
//#include "searcher.h"
//#include "data_util.h"
using namespace std;
/////////////////////////////////////////////////////////////
// gps coordinate
//
// illustrates serialization for a simple type
//
class gps_position
{
private:
	friend class boost::serialization::access;
	// When the class Archive corresponds to an output archive, the
	// & operator is defined similar to <<.  Likewise, when the class Archive
	// is a type of input archive the & operator is defined similar to >>.
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & degrees;
		ar & minutes;
		ar & seconds;
	}
	int degrees;
	int minutes;
	float seconds;
public:
	gps_position(){};
	gps_position(int d, int m, float s) :
		degrees(d), minutes(m), seconds(s)
	{}
};


struct RerankADC8 {
	int pid;
	int quantizations[8];
	template<class Archive>
	void serialize(Archive& arc, unsigned int version) {
		arc & pid;
		arc & quantizations;
	}
};

template<class T>
struct Multitable {
	/**
	*  This constructor gets width of table for each dimension
	*  @param dimensions array of sizes of table along each dimension
	*/
	//Multitable(const vector<int>& dimensions = vector<int>());
	/**
	*  This function resize the table to new dimensions
	*  @param dimensions array of sizes of table along each dimension
	*/
	//void Resize(const vector<int>& dimensions, T value = T());
	/**
	*  This function sets value in one cell
	*  @param value value to set
	*  @param cell_indices coordinates of cell in the table
	*/
	//void SetValue(T value, const vector<int>& cell_indices);
	/**
	*  This function gets value of one cell
	*  @param cell_indices coordinates of cell in the table
	*/
	//T GetValue(const vector<int>& cell_indices);
	/**
	*  Actual data as one-dimensional array
	*/
	vector<T> table;
	/**
	*  Dimensions of table
	*/
	vector<int> dimensions;
	/**
	*  Function for Boost.Serialization
	*/
	template<class Archive>
	void serialize(Archive& arc, unsigned int version) {
		arc & table;
		arc & dimensions;
	}
	/**
	*  Function converts cell coordinates to global index in a long array
	*  @param cell_indices coordinates of cell in the table
	*/
	//int GetCellGlobalIndex(const vector<int>& cell_indices) const;
};

template<class Record>
struct MultiIndex {
	vector<Record> multiindex;
	Multitable<int> cell_edges;    ///< Table with index cell edges in array
};


int main() {
	//// create and open a character archive for output
	//std::ofstream ofs("filename", std::ios::binary);
	//// create class instance
	//const gps_position g(35, 59, 24.567f);
	//// save data to archive
	//{
	//	boost::archive::binary_oarchive oa(ofs);
	//	// write class instance to archive
	//	oa << g;
	//	// archive and stream closed when destructors are called
	//}

	// ... some time later restore the class instance to its orginal state
	//gps_position newg;
	//{
	//	// create and open an archive for input
	//	std::string edge_path = "I:/Hashing/Code/OPQ_Mindex_fyt/index/fvImgDb/_cell_edges.bin";
	//	std::ifstream ifs("filename", std::ios::binary);
	//	boost::archive::binary_iarchive ia(ifs);
	//	// read class state from archive
	//	ia >> newg;
	//	// archive and stream closed when destructors are called
	//}
	
	//MultiIndexer<RerankADC8> indexer(multiplicity);
	MultiIndex<RerankADC8> multiindex_;
	std::string edge_path = "I:/Hashing/Code/OPQ_Mindex_fyt/index/fvImgDb/_cell_edges.bin";
	std::ifstream ifs(edge_path, std::ios::binary);
	boost::archive::binary_iarchive ia(ifs);
	// read class state from archive
	//gps_position newg;

	ia >> multiindex_.cell_edges;
	return 0;
}