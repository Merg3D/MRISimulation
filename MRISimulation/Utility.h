#pragma once

#include <random>
#include "Vec3.h"

#define PI 3.14159265359f
#define DEG_RAD (PI / 180.0f)

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

	}
	float get_random()
	{
		return static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX);
	}
	vec3 get_random_vec3()
	{
		float theta = get_random() * PI * 2.0f;
		float phi = get_random() * PI * 2.0f;

		float x = std::sin(theta) * std::cos(phi);
		float y = std::sin(theta) * std::sin(phi);
		float z = std::cos(theta);

		return vec3(x, y, z);
	}
};

