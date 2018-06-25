#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace boost::multiprecision;

#include <iostream>
#include <sstream>
#include <map>
using namespace std;

#include <opencv2/opencv.hpp>
using namespace cv;
#pragma comment(lib, "opencv_world340.lib")

#include "main.h"


int main(void)
{
	Mat frame = imread("picture.jpg", CV_LOAD_IMAGE_GRAYSCALE);

	if (frame.empty())
	{
		cout << "Could not read input file" << endl;
		return 1;
	}

	cout << "Allocating RAM" << endl;

	int num_wide = frame.cols;
	int num_tall = frame.rows;
	int num_dims = 3; // xyz

	vector<vector<vector<float> > > control_points;

	control_points.resize(num_wide);

	for (size_t i = 0; i < num_wide; i++)
	{
		control_points[i].resize(num_tall);

		for (size_t j = 0; j < num_tall; j++)
			control_points[i][j].resize(num_dims);
	}

	for (int j = 0; j < num_tall; j++)
	{
		for (int i = 0; i < num_wide; i++)
		{
			control_points[i][j][0] = double(i) / double(num_tall);
			control_points[i][j][1] = double(j) / double(num_tall);
			control_points[i][j][2] = double(frame.at<unsigned char>(j, i)) / 255.0;
		}
	}

	int out_num_wide = 64;// num_wide * 2;
	int out_num_tall = 48;// num_tall * 2;

	vector<int> nums;

	nums.push_back(out_num_wide);
	nums.push_back(out_num_tall);
	nums.push_back(num_wide);
	nums.push_back(num_tall);

	sort(nums.begin(), nums.end());

	int largest = nums[nums.size() - 1];

	init_fact_lut(largest);

	Mat out_frame(out_num_tall, out_num_wide, CV_8UC1);

	for (int j = 0; j < out_num_tall; j++)
	{
		for (int i = 0; i < out_num_wide; i++)
		{
			float u = i / float(out_num_wide - 1);
			float v = j / float(out_num_tall - 1);

			cout << u << " " << v << endl;

			vertex_3 v0 = bezier(u, v, control_points, num_wide, num_tall);

			out_frame.at<unsigned char>(j, i) = v0.z*255.0;
		}
	}

	imwrite("out.png", out_frame);

	return 0;
}