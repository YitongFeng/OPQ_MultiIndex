#ifndef SIFT_FV_DLL_H
#define SIFT_FV_DLL_H

#include <string>
#include <vector>
#include<opencv2/opencv.hpp>
//using namespace std;
using cv::Mat;


#ifdef DLL_EXPORT
#define DLL_OUT __declspec(dllexport)
#else 
#define DLL_OUT __declspec(dllimport)
#endif

DLL_OUT bool init_fv();
DLL_OUT bool compute_fv(Mat gray_image, Mat edge_image, Mat& fv);

#endif