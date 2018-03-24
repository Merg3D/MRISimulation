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
	void update_position(Spatial* p_spatial);
	
	// update proton (position and magnetization)
	void update_proton(Proton& p_proton);
	void update_proton_magnetization(Proton& p_proton);

	void apply_exc_pulse();
	void apply_ors_pulse();

	float B_tot(Proton& p_particle);
	float B_dip(const vec3& p_r);
	
	std::vector<Proton> protons;
	std::vector<Particle> particles;

	vec3 space_size;

	float D;
	float step_size;
	float normal_dt;
	float dt;

	float R1, R2;
	float Ms;
	float M0;
	float B0;
	float B_eq;
	float T;
	float r2;

	float particle_radius;

	float exc_pulse_flipangle; // degrees
	float ors_pulse_flipangle; // degrees

	bool running;

	float t;
	int iteration;
	int max_iterations;

	int N_protons;
	int N_particles;

	float ORS_offset;			// Hz
	float ORS_bandwidth;		// Hz

	// proton gyromagnetic ratio
	float gamma;

	Random random;

	vec3 global_magnetization;
};

