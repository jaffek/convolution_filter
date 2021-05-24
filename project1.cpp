#include"project1.h"
#include"interface.h"

int laplace_kernel, gaussian_kernel, roi_x1, roi_x2, roi_y1, roi_y2, start_calculation, gaussian_sigma, R_channel, G_channel, B_channel, choose_alg;
//laplace_kernel					define kernel size of the laplacian mask: length of the side of the matrix = 2 * laplace_kernel + 1
//gaussian_kernel					define kernel size of the gaussian mask: length of the side of the matrix = 2 * gaussian_kernel + 1
//roi_x1, roi_y1, roi_x2, roi_y2  	coordinates of the left and right bottom corner of the region of interest
//start_calculation					starts and stops the program
//gaussian_sigma				    sigma coefficient in the gaussian formula
//R_channel							define on which channel use algorythm
//G_channel							define on which channel use algorythm
//B_channel							define on which channel use algorythm
//choose_alg						define which algorythm (LoG or Gauss + Laplace) will start

Mat gauss_mask, gaussian_img, log_img, laplace_mask, ROI_img, float_img, loaded_img;
//gauss_mask						gaussian matrix
//laplace_mask						laplacian matrix
//loaded_img						the original photo loaded from the file
//float_img							loaded_img after convertion from uchar to float
//ROI_img							region of interest of the float image
//log_img							image after laplacian of gaussian operation

//declared as global variables due to the gui environment requirements in opencv (trackbars)

//***********************************************************************************************************************************
// defining initial values
void define_variables(const Mat & source_img)
{
	laplace_kernel = 0;
	gaussian_kernel = 1;
	roi_x1 = 0;
	roi_y1 = 0;
	roi_x2 = source_img.cols - 1;
	roi_y2 = source_img.rows - 1;
	start_calculation = 0;
	gaussian_sigma = 50;
	R_channel = 1;
	G_channel = 1;
	B_channel = 1;
	choose_alg = 0;
}
//------------------------------------------------------------------------------------------------------------------------------------
// calculate gaussian matrix
void create_gaussian_mask(Mat & out, int sigma_tmp, int kernel_size)
{
	float sigma = (float)sigma_tmp / 10;
	out = Mat(2 * kernel_size + 1, 2 * kernel_size + 1, CV_32FC1);
	float sum = 0, normalize = 0;
	for (int i = -kernel_size; i <= kernel_size; i++)
		for (int j = -kernel_size; j <= kernel_size; j++)
		{
			out.at<float>(i + kernel_size, j + kernel_size) = 1.0 / (2 * CV_PI*sigma*sigma)*exp(-((i * i + j * j) / (2 * sigma * sigma)));  // gaussian formula
			sum = sum + out.at<float>(i + kernel_size, j + kernel_size);
		}
	// normalization - the sum of all elements of the gaussian matrix must be equal to 1
	normalize = 1.0 / sum;                 
	for (int i = -kernel_size; i <= kernel_size; i++)
		for (int j = -kernel_size; j <= kernel_size; j++)
			out.at<float>(i + kernel_size, j + kernel_size) = out.at<float>(i + kernel_size, j + kernel_size) * normalize;
}
//------------------------------------------------------------------------------------------------------------------------------------
// calculate laplace matrix - the sum of all elements of the gaussian matrix must be equal to 0 - evereywhere 1 except the central element
void laplace(Mat & out, int kernel_size)
{
	out = Mat(2 * kernel_size + 1, 2 * kernel_size + 1, CV_32FC1);
	for (int i = -kernel_size; i <= kernel_size; i++)
		for (int j = -kernel_size; j <= kernel_size; j++)
			out.at<float>(i + kernel_size, j + kernel_size) = 1;
	out.at<float>(kernel_size, kernel_size) = -((kernel_size * 2 + 1)*(kernel_size * 2 + 1) - 1);
}
//------------------------------------------------------------------------------------------------------------------------------------
// convertion from CV_8U -> CV_32F taking into account the region of interest
void rewrite_matrix_to_float(const Mat & source_img, Mat & out, int x1, int y1, int x2, int y2, int channels_num) 
{
	if(channels_num == 3)
		out = Mat(y2 - y1, x2 - x1, CV_32FC3);
	else
		out = Mat(y2 - y1, x2 - x1, CV_32FC1);
	for (int i = x1; i < x2; i++)
		for (int j = y1; j < y2; j++)
			if (channels_num == 3)
				for (int k = 0; k < 3; k++)
					if (source_img.depth() == CV_32F)
						out.at<Vec3f>(j - y1, i - x1)[k] = source_img.at<Vec3f>(j, i)[k];
					else
						out.at<Vec3f>(j - y1, i - x1)[k] = 1.0 / 255 * source_img.at<Vec3b>(j, i)[k];
			else
				if (source_img.depth() == CV_32F)
					out.at<float>(j - y1, i - x1) = source_img.at<float>(j, i);
				else
					out.at<float>(j - y1, i - x1) = 1.0 / 255 * source_img.at<uchar>(j, i);
}
//------------------------------------------------------------------------------------------------------------------------------------
// convolution of source image and mask
void convolution_img(const Mat & source_img, Mat & out, Mat mask)
{
	int kernel_size = (mask.rows - 1) / 2;
	if (source_img.channels() == 1)
	{
		out = Mat(source_img.rows - 2 * kernel_size, source_img.cols - 2 * kernel_size, CV_32FC1);
		for (int i = kernel_size; i < source_img.rows - kernel_size; i++)
			for (int j = kernel_size; j < source_img.cols - kernel_size; j++)                     // navigating the source image
			{
				float sum = 0;
				for (int k = -kernel_size; k <= kernel_size; k++)
					for (int m = -kernel_size; m <= kernel_size; m++)                             // navigating the kernel
					{
						float component = 0;
						component = source_img.at<float>(i + k, j + m) * mask.at<float>(-k + kernel_size, -m + kernel_size);
						sum += component;
					}
				out.at<float>(i - kernel_size, j - kernel_size) = sum;
			}
	}
	else
	{
		out = Mat(source_img.rows - 2 * kernel_size, source_img.cols - 2 * kernel_size, CV_32FC3);
		for (int i = kernel_size; i < source_img.rows - kernel_size; i++)
			for (int j = kernel_size; j < source_img.cols - kernel_size; j++)                   // navigating the source image
			{
				float sum[3] = { 0,0,0 };
				for (int k = -kernel_size; k <= kernel_size; k++)
					for (int m = -kernel_size; m <= kernel_size; m++)                           // navigating the kernel
					{
						float component[3] = { 0,0,0 };
						for (int n = 0; n < 3; n++)
						{
							component[n] = source_img.at<Vec3f>(i + k, j + m)[n] * mask.at<float>(-k + kernel_size, -m + kernel_size);
							sum[n] += component[n];
						}
					}
				if(R_channel == 1)
					out.at<Vec3f>(i - kernel_size, j - kernel_size)[2] = sum[2];
				if (G_channel == 1)
					out.at<Vec3f>(i - kernel_size, j - kernel_size)[1] = sum[1];
				if (B_channel == 1)
					out.at<Vec3f>(i - kernel_size, j - kernel_size)[0] = sum[0];
			}
	}
}
//------------------------------------------------------------------------------------------------------------------------------------
// another implementation - calculation of tha Laplacian of Gaussian mask using known LoG formula
void laplacian_of_gaussian(Mat & out, int sigma_tmp, int kernel_size)
{
	float sigma = (float)sigma_tmp / 10;
	float max_element = 0;
	out = Mat(2 * kernel_size + 1, 2 * kernel_size + 1, CV_32FC1);
	float sum = 0, normalize_mask = 0;
	for (int i = -kernel_size; i <= kernel_size; i++)
		for (int j = -kernel_size; j <= kernel_size; j++)
		{
			out.at<float>(i + kernel_size, j + kernel_size) = -1.0 / (CV_PI*pow(sigma,4))*(float)(1-(float)(i*i+j*j)/(2*sigma*sigma))*exp(-((float)(i * i + j * j) / (2 * sigma * sigma)));  // laplacian of gaussian formula
			sum += out.at<float>(i + kernel_size, j + kernel_size);
		}
	//  normalize - sum have to be 0
	normalize_mask = -sum / (pow((2 * kernel_size + 1), 2));
	for (int i = -kernel_size; i <= kernel_size; i++)
		for (int j = -kernel_size; j <= kernel_size; j++)
		{
			out.at<float>(i + kernel_size, j + kernel_size) = out.at<float>(i + kernel_size, j + kernel_size) + normalize_mask;
			if (abs(out.at<float>(i + kernel_size, j + kernel_size)) > max_element)
				max_element = abs(out.at<float>(i + kernel_size, j + kernel_size));
		}
	sum = 0;
	for (int i = -kernel_size; i <= kernel_size; i++)
		for (int j = -kernel_size; j <= kernel_size; j++)
		{
			out.at<float>(i + kernel_size, j + kernel_size) = out.at<float>(i + kernel_size, j + kernel_size) / max_element;
			sum += out.at<float>(i + kernel_size, j + kernel_size);
		}
}
//------------------------------------------------------------------------------------------------------------------------------------
// function called when user move the START/STOP trackbar
void on_trackbar(int, void*)
{
	if (start_calculation == 1)
	{
		if (roi_x1 < roi_x2 && roi_y1 < roi_y2)
			rewrite_matrix_to_float(float_img, ROI_img, roi_x1, roi_y1, roi_x2, roi_y2,float_img.channels());
		if (choose_alg == 1)
		{
			laplacian_of_gaussian(laplace_mask, gaussian_sigma, laplace_kernel + 1);    // use with laplacian_of_gaussian method
			convolution_img(ROI_img, log_img, laplace_mask);							// use with laplacian_of_gaussian method
		}
		else
		{
			create_gaussian_mask(gauss_mask, gaussian_sigma, gaussian_kernel);  // gaussian + laplacian method
			laplace(laplace_mask, laplace_kernel+1);							// gaussian + laplacian method
			convolution_img(ROI_img, gaussian_img, gauss_mask);				    // gaussian + laplacian method
			convolution_img(gaussian_img, log_img, laplace_mask);				// gaussian + laplacian method
		}
		display(1, log_img);
	}
	else
		display(0, log_img);
}
//------------------------------------------------------------------------------------------------------------------------------------
// calls other functions instead main function in main.cpp
void main_function()
{
	loaded_img = imread("lena1.jpg");             // RGB image
	Mat trackbars_image = imread("trackbars.jpg");
	//loaded_img = imread("cheetah.jpg");				// grayscale image
	rewrite_matrix_to_float(loaded_img, float_img, 0, 0, loaded_img.cols, loaded_img.rows,loaded_img.channels());
	define_variables(float_img);
	trackbars(loaded_img.rows, loaded_img.cols, float_img,trackbars_image);
	waitKey(0);
}
//***********************************************************************************************************************************