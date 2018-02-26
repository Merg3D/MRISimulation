#pragma once

#include <vector>

#include "Vec3.h"

struct Particle
{
	vec3 position;
};

struct SPIO : public Particle
{

};

struct Proton : public Particle
{
	Proton()
	{
		magnetization = vec3(0.0f, 0.0f, 1.0f);
	}

	vec3 magnetization;
};

class Simulator
{
public:
	Simulator();
	~Simulator();

	void iterate();

	float get_signal();
	float get_contrast();

private:
	// calculate global magnetization
	void calculate_global_magnetization();

	// assert periodic boundary conditions
	void assert_in_space(Particle* p_particle);

	// update particle (position and magnetization) for next iteration
	void update_particle(Particle* p_particle);
	void update_particle_position(Particle* p_particle);
	void update_particle_magnetization(Particle* p_particle);

	std::vector<Particle*> particles;

	vec3 space_size;

	float D;
	float step_size;
	float normal_dt;
	float dt;

	float R1, R2;
	float T1;

	vec3 global_magnetization;
};

