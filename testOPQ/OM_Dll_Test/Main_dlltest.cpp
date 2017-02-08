#include <iostream>
#include <vector>
#include <string>
#include "OPQ_MI_DLL.h"
using namespace std;

int main(){
	init();
	vector<pair<string, float>> query_result;
	string path = "T:/20161108/16/total/P1595.jpg";
	query(path.c_str(), query_result);

	return 0;
}