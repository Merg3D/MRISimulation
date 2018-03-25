#include "Simulator.h"

#include <fstream>

#include "Mat4.h"


Simulator::Simulator()
{
	max_iterations = 20;

	double Cc = 1.0e-3;						// M

	N_protons = 5e5;
	N_particles = 1e2;

	D = 3.0e-9;							// m^2 s^-1

	double T1 = 3.0;						// s
	double T2 = 80.0e-3;					// s

	R1 = 1.0 / T1;							// Hz
	R2 = 1.0 / T2;							// Hz

	Ms = 100.0;							// Am^2 / kg [Fe]
	M0 = 1.0;

	double mu_0 = 4e-7 * PI;
	B_eq = 0.72 * 5185 * mu_0 / 3.0 * Ms; // Tesla (0.156)

	T = 310.0;								// K

	particle_radius = 20e-9;				// m

	gamma = 267.513e6;						// rad s^−1 T^−1)

	normal_dt = pow(particle_radius, 2) / (6 * D);	// s
	dt = normal_dt;							// s

	exc_pulse_flipangle = 10.0;			// degrees
	ors_pulse_flipangle = 90.0;			// degrees

	dt = 0.1;

	space_size = vec3d(1.5e-5);				// m

	B0 = 7.0;								// T

	ORS_frequency = 0.0;
	ORS_bandwidth = 5000.0;
}

Simulator::~Simulator()
{
	protons.clear();
	particles.clear();
}

void Simulator::start()
{
	init();
	apply_exc_pulse();

	apply_ors_pulse();

	running = true;

	while (running)
	{
		if (iteration >= max_iterations)
		{
			running = false;
			break;
		}

		double signal = get_signal();
		signals.push_back(signal);
		std::cout << "Iteration: " << iteration << ", time: " << t << ", signal = " << signal << std::endl;

		iterate();

		t += dt;
		iteration++;
	}

	output();
}

void Simulator::init()
{
	protons.resize(N_protons);
	particles.resize(N_particles);

	for (int c = 0; c < N_protons; c++)
	{
		protons[c].position = vec3d(random.get_random(), random.get_random(), random.get_random()) * space_size;
		protons[c].magnetization = vec3d(0, 0, M0);
	}

	for (int c = 0; c < N_particles; c++)
	{
		particles[c].position = vec3d(random.get_random(), random.get_random(), random.get_random()) * space_size;
	}
}

void Simulator::iterate()
{
	step_size = std::sqrt(6 * D * dt);

	for (int c = 0; c < N_protons; c++)
	{
		Proton& proton = protons[c];

		update_proton(proton);
		assert_in_space(&proton);
	}
}

double Simulator::get_signal()
{
	calculate_global_magnetization();
	vec3d M = global_magnetization;
	return std::sqrt(M.x * M.x + M.y * M.y);
}

double Simulator::get_contrast()
{
	return 0.0;
}

void Simulator::calculate_global_magnetization()
{
	vec3d sum = 0.0;

	for (int c = 0; c < N_protons; c++)
		sum += protons[c].magnetization;

	global_magnetization = sum / (double) N_protons;
}

void Simulator::assert_in_space(Spatial* p_particle)
{
	if (p_particle->position.x > space_size.x)
		p_particle->position.x -= space_size.x * 2.0;
	else if (p_particle->position.x < space_size.x)
		p_particle->position.x += space_size.x * 2.0;

	if (p_particle->position.y > space_size.y)
		p_particle->position.y -= space_size.y * 2.0;
	else if (p_particle->position.y < space_size.y)
		p_particle->position.y += space_size.y * 2.0;

	if (p_particle->position.z > space_size.z)
		p_particle->position.z -= space_size.z * 2.0;
	else if (p_particle->position.z < space_size.z)
		p_particle->position.z += space_size.z * 2.0;
}

void Simulator::update_proton(Proton& p_particle)
{
	update_position(&p_particle);
	update_proton_magnetization(p_particle);
}

void Simulator::update_position(Spatial* p_particle)
{
	p_particle->position += random.get_random_vec3() * step_size;
}

void Simulator::update_proton_magnetization(Proton& p_proton)
{
	vec3d& magn = p_proton.magnetization;

	double B_t = B_tot(p_proton);
	double theta = gamma * B_t * dt;
	double E1 = std::exp(-R1 * dt);
	double E2 = std::exp(-R2 * dt);

	magn.x = E1 * (magn.x * cos(theta) + magn.y * sin(theta));
	magn.y = E1 * (magn.x * sin(theta) - magn.y * cos(theta));
	magn.z = E2 * magn.z + (1.0 - std::exp(-dt * R1)) * M0;
}

void Simulator::apply_exc_pulse()
{
	mat4 rot;
	rot.rotate_x(exc_pulse_flipangle * DEG_RAD);

	for (int c = 0; c < N_protons; c++)
		protons[c].magnetization = (rot * vec4d(protons[c].magnetization, 1.0)).get_xyz();
}

void Simulator::apply_ors_pulse()
{
	mat4 rot;
	rot.rotate_x(ors_pulse_flipangle * DEG_RAD);

	std::cout << "Applying ORS pulse" << std::endl;

	int count = 0;

	for (int c = 0; c < N_protons; c++)
	{
		vec3d& magn = protons[c].magnetization;

		double B_t = B_tot(protons[c]);
		double offset = gamma * B_t - gamma * B0; // relative larmor frequency

		if (offset > ORS_frequency - ORS_bandwidth / 2.0 && offset < ORS_frequency + ORS_bandwidth / 2.0)
		{
			protons[c].magnetization = (rot * vec4d(protons[c].magnetization, 1.0)).get_xyz();
			count++;
		}

		// store frequency
		offsets.push_back(offset);
	}

	std::cout << "Finished ORS pulse, saturated " << (double) count / N_protons * 100 << "% of the protons" << std::endl;
}

double Simulator::B_tot(Proton& p_proton)
{
	double sum = 0.0;

	for (int i = 0; i < N_particles; i++)
		sum += B_dip(p_proton.position - particles[i].position);

	return B0 + sum;
}

double Simulator::B_dip(const vec3d& p_r)
{
	double theta = vec3d(0, 0, 1).angle(p_r);
	return B_eq * pow(particle_radius / p_r.length(), 3) * (3 * pow(cos(theta), 2) - 1.0);
}

void Simulator::output()
{
	std::ofstream file_off;
	std::ofstream file_sig;
	file_off.open("offsets.txt");
	file_sig.open("signals.txt");

	for (int c = 0; c < offsets.size(); c++)
		file_off << offsets[c] << "\n";

	for (int c = 0; c < signals.size(); c++)
		file_sig << signals[c] << "\n";
	
	file_off.close();
	file_sig.close();
}
