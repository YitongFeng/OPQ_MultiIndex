#ifndef OPQ_MI_DLL_H
#define OPQ_MI_DLL_H

#include <string>
#include <vector>
using namespace std;

#define OPQ_MI_DLL __declspec(dllexport)
OPQ_MI_DLL bool init();
OPQ_MI_DLL bool release();
//DLL_OUT bool query_fv(char * filepath);
OPQ_MI_DLL bool query(char * filepath, vector<pair<string, float>>& query_result, bool proposal_flag, vector<int> class_choose);

#endif

