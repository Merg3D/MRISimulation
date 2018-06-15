#pragma once

#include <vector>
#include <thread>
#include <atomic>

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

// Experiment: class describing all parameters for the experiment + contains the results
struct Experiment
{
	double dt;

	int max_iterations = 1000;

	int N_protons = 1e4;
	int N_particles = 1e2;

	double ors_frequency = 200.0;
	double ors_bandwidth = 100.0;

	double Cc;

	vec3d volume = vec3d(1e-6); // m

	int id = -1;

	int averages = 4;

	std::vector<double> results;
};

class Simulator
{
public:
	Simulator();
	~Simulator();

	void start();

	void start_experiment(Experiment& p_exp);

	// calculate and store signal
	double get_signal();

private:

	void init();
	void iterate();

	// not really necessary, but helpful for debugging
	double get_z_magnetization();

	// calculate global magnetization
	void calculate_global_magnetization();

	// assert periodic boundary conditions
	void assert_in_space(Proton& p_proton);
	
	// update position for next iteration
	void update_proton_positions(int start, int end);
	void update_proton_position(Proton& p_proton);
	
	// update proton (position and magnetization)
	void update_proton_magnetizations(int start, int end);
	void update_proton_magnetization(Proton& p_proton);
	
	// apply excitation and off-resonance pulses along x-axis
	void apply_exc_pulse();
	void apply_ors_pulse();

	// calculate total magnetic field
	double B_tot(Proton& p_particle);

	// calculate dipole field for SPIO particle given distance to proton 
	double B_dip(const vec3d& p_r);

	// write results to text files
	void save();
	
	std::vector<Proton> protons;
	std::vector<Particle> particles;

	// results
	std::vector<double> signals;
	std::vector<double> offsets;

	std::vector<Experiment> experiments;
	std::vector<std::thread> threads;

	int current_exp;

	vec3d volume;

	double D;
	double step_size;
	double normal_dt;
	double dt;

	double R1, R2;
	double Ms;
	double M0;
	double B0;
	double B1;
	double B_eq;
	double T;					// temperature

	double particle_radius;

	double exc_pulse_flipangle; // degrees
	double ors_pulse_flipangle; // degrees

	bool running;
	bool applying_gradient;

	double t;
	int iteration;
	int max_iterations;

	int N_protons;
	int N_particles;

	int N_threads;

	double ors_frequency;		// Hz
	double ors_bandwidth;		// Hz

	double exc_pulse_time;
	double gradient_duration;

	bool has_applied_exc_pulse;

	// proton gyromagnetic ratio
	double gamma;

	Random random;

	vec3d global_magnetization;
};

