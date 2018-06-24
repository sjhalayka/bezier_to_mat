#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace boost::multiprecision;

#include <iostream>
#include <sstream>
#include <map>
using namespace std;

#include <opencv2/opencv.hpp>
using namespace cv;
#pragma comment(lib, "opencv_world340.lib")

#include "primitives.h"

double mp_to_double(cpp_dec_float_100 &b)
{
	ostringstream oss;
	oss << b;

	double ret;
	istringstream iss(oss.str());
	iss >> ret;

	return ret;
}

long long unsigned mp_to_int(cpp_dec_float_100 &b)
{
	ostringstream oss;
	oss << b;

	long long unsigned ret;
	istringstream iss(oss.str());
	iss >> ret;

	return ret;
}

vector<cpp_dec_float_100> fact_lut;

void init_fact_lut(long long unsigned int n_max)
{
	fact_lut.resize(n_max + 1);

	for (long long unsigned int i = 0; i <= n_max; i++)
	{
		cpp_dec_float_100 ret = 1;

		for (cpp_dec_float_100 k = i; k > 0; k--)
			ret *= k;

		fact_lut[i] = ret;
	}
}

cpp_dec_float_100 fact(cpp_dec_float_100 n)
{
	return fact_lut[mp_to_int(n)];
}

map<pair<cpp_dec_float_100, cpp_dec_float_100>, cpp_dec_float_100> binomial_cache;

cpp_dec_float_100 binomial(cpp_dec_float_100 n, cpp_dec_float_100 k)
{
	pair<cpp_dec_float_100, cpp_dec_float_100> p(n, k);

	map<pair<cpp_dec_float_100, cpp_dec_float_100>, cpp_dec_float_100>::const_iterator ci = binomial_cache.find(p);

	if (ci != binomial_cache.end())
		return ci->second;

	//     n!
	// -----------
	// k! (n - k)!

	cpp_dec_float_100 ret = fact(n) / (fact(k)*fact(n - k));

	binomial_cache[p] = ret;

	return ret;
}

map<pair<cpp_dec_float_100, cpp_dec_float_100>, cpp_dec_float_100> pow_cache;

cpp_dec_float_100 cached_pow(cpp_dec_float_100 &base, cpp_dec_float_100 &exponent)
{
	pair<cpp_dec_float_100, cpp_dec_float_100> p(base, exponent);

	map<pair<cpp_dec_float_100, cpp_dec_float_100>, cpp_dec_float_100>::const_iterator ci = pow_cache.find(p);

	if (ci != pow_cache.end())
		return ci->second;

	cpp_dec_float_100 ret = pow(base, exponent);

	pow_cache[p] = ret;

	return ret;
}


vertex_3 bezier(const double u, const double v, vector<vector<vector<float> > > &control_points, size_t num_wide, size_t num_tall)
{
	vertex_3 ret;

	cpp_dec_float_100 ret_x;
	cpp_dec_float_100 ret_y;
	cpp_dec_float_100 ret_z;

	cpp_dec_float_100 u_long = u;
	cpp_dec_float_100 v_long = v;

	const size_t Ni = num_wide - 1;
	const size_t Nj = num_tall - 1;

	for (size_t i = 0; i <= Ni; i++)
	{
		for (size_t j = 0; j <= Nj; j++)
		{
			cpp_dec_float_100 binpow_u = binomial(Ni, i) * cached_pow(u_long, cpp_dec_float_100(i)) * cached_pow(cpp_dec_float_100(1 - u_long), cpp_dec_float_100(Ni - i));
			cpp_dec_float_100 binpow_v = binomial(Nj, j) * cached_pow(v_long, cpp_dec_float_100(j)) * cached_pow(cpp_dec_float_100(1 - v_long), cpp_dec_float_100(Nj - j));

			ret_x += binpow_u * binpow_v * control_points[i][j][0];
			ret_y += binpow_u * binpow_v * control_points[i][j][1];
			ret_z += binpow_u * binpow_v * control_points[i][j][2];
		}
	}

	ret.x = mp_to_double(ret_x);
	ret.y = mp_to_double(ret_y);
	ret.z = mp_to_double(ret_z);

	return ret;
}



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

	int out_num_wide = 16;// num_wide * 2;
	int out_num_tall = 12;// num_tall * 2;

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