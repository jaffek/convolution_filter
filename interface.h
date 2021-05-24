#ifndef INTERFACE_H 
#define INTERFACE_H
#include<opencv2/opencv.hpp>
#include<iostream>
#include "project1.h"

using namespace cv;

//------------------------------------------------------------------------------------------------------------------------------------
// display window with effect or destroy
void display(int open_close_window, const Mat & image_to_display);
//------------------------------------------------------------------------------------------------------------------------------------
// create trackbars to edit settings of the Laplacian of Gaussian transformation
void trackbars(int img_rows, int img_cols, Mat source_img, Mat trackbars_img);

#endif
