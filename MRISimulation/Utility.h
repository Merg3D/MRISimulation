#pragma once

#include <random>
#include <time.h>

#include "Vec3.h"

#define PI 3.14159265359f
#define DEG_RAD (PI / 180.0)

class Utility
{
public:
	Utility();
	~Utility();
};

class Random
{
public:
	Random()
	{
		/* initialize random seed: */
		srand(time(NULL));
	}
	double get_random()
	{
		return static_cast<double>(std::rand()) / static_cast<double>(RAND_MAX);
	}
	vec3d get_random_vec3()
	{
		double theta = get_random() * PI * 2.0;
		double phi = get_random() * PI * 2.0;

		double x = std::sin(theta) * std::cos(phi);
		double y = std::sin(theta) * std::sin(phi);
		double z = std::cos(theta);

		return vec3d(x, y, z);
	}
};

