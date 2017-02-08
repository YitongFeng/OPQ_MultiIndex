#ifndef _H_CNN_API_H
#define _H_CNN_API_H

#include <opencv2/opencv.hpp>

#define DLL_EXPORT	__declspec(dllexport)


using cv::Mat;

DLL_EXPORT void cnn_init();
DLL_EXPORT void cnn_fea(Mat in, Mat& out);

#endif