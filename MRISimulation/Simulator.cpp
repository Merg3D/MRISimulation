#include "Simulator.h"



Simulator::Simulator()
{
	max_iterations = 1e3;

	N_protons = 1e8;
	N_particles = 1e4;

	D = 3.0e-9f;				// m^2 s^-1

	float T1 = 3.0f;			// s
	float T2 = 80.0e-3f;		// s

	R1 = 1.0f / T1;				// Hz
	R2 = 1.0f / T2;				// Hz

	Mz = 100.0f;				// Am^2 / kg [Fe]
	B_eq = 16.0f;				// T
	T = 300.0f;					// K
	r2 = 600.0f;				// mM^-1 s^-1

	particle_radius = 20e-9f;	// m
	proton_radius = 0.87e-9f;	// m


}


Simulator::~Simulator()
{
}

void Simulator::start()
{
	init();

	running = true;

	while (running)
	{
		if (t >= max_iterations)
			running = false;

		iterate();

		t++;
	}
}

void Simulator::init()
{
	protons.reserve(N_protons);
	particles.reserve(N_particles);

	for (int c = 0; c < N_protons; c++)
	{
		protons[c].position = vec3(random.get_random(), random.get_random(), random.get_random()) * space_size;
		protons[c].magnetization = vec3(0, 0, Mz);
	}

	for (int c = 0; c < N_particles; c++)
	{
		particles[c].position = vec3(random.get_random(), random.get_random(), random.get_random());
	}
}

void Simulator::iterate()
{
	step_size = std::sqrt(6 * D * dt);

	for (int c = 0; c < particles.size(); c++)
	{
		Particle& p = particles[c];

		update_particle(p);
		assert_in_space(&p);

		std::cout << "Iteration: " << c << std::endl;

	}
}

float Simulator::get_signal()
{
	vec3 M = global_magnetization;
	return std::sqrt(M.x * M.x + M.y * M.y);
}

float Simulator::get_contrast()
{
	return 0.0f;
}

void Simulator::calculate_global_magnetization()
{
	vec3 sum = 0.0f;

	for (int c = 0; c < particles.size(); c++)
	{
		Proton& p = protons[c];
		sum += p.magnetization;
	}

	global_magnetization = sum / (float)particles.size();
}

void Simulator::assert_in_space(Spatial* p_particle)
{
	if (p_particle->position.x > space_size.x)
		p_particle->position.x -= space_size.x * 2.0f;
	else if (p_particle->position.x < space_size.x)
		p_particle->position.x += space_size.x * 2.0f;

	if (p_particle->position.y > space_size.y)
		p_particle->position.y -= space_size.y * 2.0f;
	else if (p_particle->position.y < space_size.y)
		p_particle->position.y += space_size.y * 2.0f;

	if (p_particle->position.z > space_size.z)
		p_particle->position.z -= space_size.z * 2.0f;
	else if (p_particle->position.z < space_size.z)
		p_particle->position.z += space_size.z * 2.0f;
}

void Simulator::update_particle(Particle& p_particle)
{
	update_position(&p_particle);
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

void Simulator::update_proton_magnetization(Proton& p_particle)
{
	vec3& magn = p_particle.magnetization;

	float gamma = 0.0f;
	float B_tot = 0.0f;
	float theta = gamma * B_tot * (t + dt) * dt;
	float E1 = std::exp(-R2 * dt);
	float E2 = std::exp(-R1 * dt);

	magn.x += E1 * (magn.x * std::cos(theta) + magn.y * std::sin(theta));
	magn.y += E1 * (magn.x * std::sin(theta) - magn.y * std::cos(theta));
	magn.z += E2 * magn.z + (1.0f - std::exp(-dt * R1));
}
