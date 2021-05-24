#include "interface.h"

//------------------------------------------------------------------------------------------------------------------------------------
// display window with effect or destroy
void display(int open_close_window, const Mat & image_to_display)
{
	if (open_close_window == 1)
	{
		namedWindow("Edges", WINDOW_NORMAL);
		imshow("Edges", image_to_display);
	}
	else
		destroyWindow("Edges");
}
//------------------------------------------------------------------------------------------------------------------------------------
// create trackbars to edit settings of the Laplacian of Gaussian transformation
void trackbars(int img_rows, int img_cols, Mat source_img, Mat trackbars_img)
{
	namedWindow("settings", CV_WINDOW_NORMAL);
	namedWindow("Original", CV_WINDOW_NORMAL);
	createTrackbar("Sigma *10", "settings", &gaussian_sigma, 100);
	createTrackbar("Gkernel", "settings", &gaussian_kernel, 10);
	createTrackbar("Lkernel +1", "settings", &laplace_kernel, 5);
	createTrackbar("R channel", "settings", &R_channel, 1);
	createTrackbar("G channel", "settings", &G_channel, 1);
	createTrackbar("B channel", "settings", &B_channel, 1);
	createTrackbar("ROI X1", "settings", &roi_x1, img_cols - 2);
	createTrackbar("ROI Y1", "settings", &roi_y1, img_rows - 2);
	createTrackbar("ROI X2", "settings", &roi_x2, img_cols - 1);
	createTrackbar("ROI Y2", "settings", &roi_y2, img_rows - 1);
	createTrackbar("G+L/LoG", "settings", &choose_alg, 1);
	createTrackbar("START/STOP", "settings", &start_calculation, 1, on_trackbar);
	imshow("settings", trackbars_img);
	imshow("Original", source_img);
}