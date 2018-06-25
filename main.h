#pragma once

#include "primitives.h"


inline double mp_to_double(cpp_dec_float_100 &b)
{
	ostringstream oss;
	oss << b;

	double ret;
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

map<pair<long long unsigned int, long long unsigned int>, cpp_dec_float_100> binomial_cache;

cpp_dec_float_100 binomial(long long unsigned int n, long long unsigned int k)
{
	pair<long long unsigned int, long long unsigned int> p(n, k);

	map<pair<long long unsigned int, long long unsigned int>, cpp_dec_float_100>::const_iterator ci = binomial_cache.find(p);

	if (ci != binomial_cache.end())
		return ci->second;

	//     n!
	// -----------
	// k! (n - k)!

	cpp_dec_float_100 ret = fact_lut[n] / (fact_lut[k]*fact_lut[n - k]);

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
