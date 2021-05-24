#ifndef PROJECT1_H 
#define PROJECT1_H
#include<opencv2/opencv.hpp>
#include<iostream>

using namespace cv;
using namespace std;

extern int laplace_kernel, gaussian_kernel, roi_x1, roi_x2, roi_y1, roi_y2, key, start_calculation, gaussian_sigma, R_channel, G_channel, B_channel, choose_alg;
//laplace_kernel					define kernel size of the laplacian mask: length of the side of the matrix = 2 * laplace_kernel + 1
//gaussian_kernel					define kernel size of the gaussian mask: length of the side of the matrix = 2 * gaussian_kernel + 1
//roi_x1, roi_y1, roi_x2, roi_y2  	coordinates of the left and right bottom corner of the region of interest
//start_calculation					starts and stops the program
//gaussian_sigma				    sigma coefficient in the gaussian formula
//R_channel							define on which channel use algorythm
//G_channel							define on which channel use algorythm
//B_channel							define on which channel use algorythm
//choose_alg						define which algorythm (LoG or Gauss + Laplace) will start

extern Mat gauss_mask, gaussian_img, log_img, laplace_mask, ROI_img, float_img, loaded_img;
//gauss_mask						gaussian matrix
//laplace_mask						laplacian matrix
//loaded_img						the original photo loaded from the file
//float_img							loaded_img after convertion from uchar to float
//ROI_img							region of interest of the float image
//log_img							image after laplacian of gaussian operation

//declared as global variables due to the gui environment requirements in opencv (trackbars)

//***********************************************************************************************************************************
// defining initial values
void define_variables();
//------------------------------------------------------------------------------------------------------------------------------------
// calculate gaussian matrix
void create_gaussian_mask(Mat & out, int sigma_tmp, int kernel_size);
//------------------------------------------------------------------------------------------------------------------------------------
// convertion from CV_8U -> CV_32F taking into account the region of interest
void rewrite_matrix_to_float(const Mat & source_img, Mat & out, int x1, int y1, int x2, int y2);
//------------------------------------------------------------------------------------------------------------------------------------
// convertion from CV_8U -> CV_32F
void uchar_to_float_conv(const Mat & source_img, Mat & conv_img);
//------------------------------------------------------------------------------------------------------------------------------------
// convolution of source image and mask
void convolution_img(const Mat & source_img, Mat & out, Mat mask);
//------------------------------------------------------------------------------------------------------------------------------------
// another implementtion - calculation of tha Laplacian of Gaussian mask using a one LoG formula
void laplacian_of_gaussian(Mat & out, int sigma_tmp, int kernel_size);
//------------------------------------------------------------------------------------------------------------------------------------
// function called when user move the START/STOP trackbar
void on_trackbar(int, void*);
//------------------------------------------------------------------------------------------------------------------------------------
// calls other functions
void main_function();
//***********************************************************************************************************************************
#endif