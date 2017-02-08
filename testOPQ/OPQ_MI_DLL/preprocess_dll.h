#ifndef PREPROCESS_DLL_H
#define PREPROCESS_DLL_H
#include<opencv2/opencv.hpp>
//using namespace std;
using cv::Mat;

#ifdef DLL_EXPORT
#define DLL_OUT __declspec(dllexport)
#else 
#define DLL_OUT __declspec(dllimport)
#endif

DLL_OUT bool init_preprocess();
DLL_OUT bool preprocess(Mat image_in, Mat& image_edge_out, Mat& image_out);

#endif