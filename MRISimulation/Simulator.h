#pragma once

#include <vector>

#include "Vec3.h"
#include "Utility.h"

struct Spatial
{
	vec3d position;
};

struct Particle : public Spatial
{

};

struct Proton : public Spatial
{
	vec3d magnetization;
};

class Simulator
{
public:
	Simulator();
	~Simulator();

	void start();

	double get_signal();
	double get_contrast();

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

	// apply excitation and off-resonance pulses
	void apply_exc_pulse();
	void apply_ors_pulse();

	// calculate total magnetic field
	double B_tot(Proton& p_particle);

	// calculate dipole field for SPIO particle given distance to proton 
	double B_dip(const vec3d& p_r);

	void output();
	
	std::vector<Proton> protons;
	std::vector<Particle> particles;

	// results
	std::vector<double> signals;
	std::vector<double> offsets;

	vec3d space_size;

	double D;
	double step_size;
	double normal_dt;
	double dt;

	double R1, R2;
	double Ms;
	double M0;
	double B0;
	double B_eq;
	double T;					// temperature

	double particle_radius;

	double exc_pulse_flipangle; // degrees
	double ors_pulse_flipangle; // degrees

	bool running;

	double t;
	int iteration;
	int max_iterations;

	int N_protons;
	int N_particles;

	double ORS_frequency;		// Hz
	double ORS_bandwidth;		// Hz

	// proton gyromagnetic ratio
	double gamma;

	Random random;

	vec3d global_magnetization;
};

