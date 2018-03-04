#pragma once

#include <vector>

#include "Vec3.h"
#include "Utility.h"

struct Spatial
{
	vec3 position;
};

struct Particle : public Spatial
{

};

struct Proton : public Spatial
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

	void start();

	float get_signal();
	float get_contrast();

private:

	void init();
	void iterate();

	// calculate global magnetization
	void calculate_global_magnetization();

	// assert periodic boundary conditions
	void assert_in_space(Spatial* p_particle);
	
	// update position for next iteration
	void update_position(Spatial* p_particle);

	// update particle
	void update_particle(Particle& p_particle);

	// update proton (position and magnetization)
	void update_proton(Proton& p_particle);
	void update_proton_magnetization(Proton& p_particle);
	
	std::vector<Proton> protons;
	std::vector<Particle> particles;

	vec3 space_size;

	float D;
	float step_size;
	float normal_dt;
	float dt;

	float R1, R2;
	float Mz;
	float B_eq;
	float T;
	float r2;
	float particle_radius;
	float proton_radius;

	bool running;

	int t;
	int max_iterations;

	int N_protons;
	int N_particles;

	Random random;

	vec3 global_magnetization;
};

