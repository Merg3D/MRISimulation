#include "Simulator.h"


#include "Mat4.h"


Simulator::Simulator()
{
	max_iterations = 1e2;

	float avogadro = 6.02214e23f;
	float Cc = 1.0e-3f;						// M

	N_protons = 1e5;
	N_particles = 1;

	D = 3.0e-9f;							// m^2 s^-1

	float T1 = 3.0f;						// s
	float T2 = 80.0e-3f;					// s

	R1 = 1.0f / T1;							// Hz
	R2 = 1.0f / T2;							// Hz

	Ms = 100.0f;							// Am^2 / kg [Fe]
	M0 = 1.0f;

	float mu_0 = 4e-7 * PI;
	B_eq = 0.72f * 5185 * mu_0 / 3.0f * Ms; // Tesla (0.156)

	T = 310.0f;								// K
	r2 = 600.0f;							// mM^-1 s^-1

	particle_radius = 20e-9f;				// m

	gamma = 267.513e6f;						// rad s^−1 T^−1)

	normal_dt = pow(particle_radius, 2) / (6 * D);	// s
	dt = normal_dt;							// s

	exc_pulse_flipangle = 10.0f;			// degrees
	ors_pulse_flipangle = 90.0f;			// degrees

	dt = 1;

	space_size = vec3(1e-3f);				// m

	B0 = 7.0f;								// T
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
			running = false;

		iterate();

		t += dt;
		iteration++;

		std::cout << "Iteration: " << iteration << ", time: " << t << ", signal = " << get_signal() << std::endl;
	}
}

void Simulator::init()
{
	protons.resize(N_protons);
	particles.resize(N_particles);

	for (int c = 0; c < N_protons; c++)
	{
		protons[c].position = vec3(random.get_random(), random.get_random(), random.get_random()) * space_size;
		protons[c].magnetization = vec3(0, 0, M0);
	}

	for (int c = 0; c < N_particles; c++)
	{
		particles[c].position = vec3(random.get_random(), random.get_random(), random.get_random()) * space_size;
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

float Simulator::get_signal()
{
	calculate_global_magnetization();
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

	for (int c = 0; c < N_protons; c++)
		sum += protons[c].magnetization;

	global_magnetization = sum / (float) N_protons;
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
	vec3& magn = p_proton.magnetization;

	float B_t = B_tot(p_proton);
	float theta = gamma * B_t * dt;
	float E1 = std::exp(-R1 * dt);
	float E2 = std::exp(-R2 * dt);

	magn.x = E1 * (magn.x * cos(theta) + magn.y * sin(theta));
	magn.y = E1 * (magn.x * sin(theta) - magn.y * cos(theta));
	magn.z = E2 * magn.z + (1.0f - std::exp(-dt * R1)) * M0;
}

void Simulator::apply_exc_pulse()
{
	mat4 rot;
	rot.rotate_x(exc_pulse_flipangle * DEG_RAD);

	for (int c = 0; c < N_protons; c++)
		protons[c].magnetization = (rot * vec4(protons[c].magnetization, 1.0f)).get_xyz();
}

void Simulator::apply_ors_pulse()
{
	mat4 rot;
	rot.rotate_x(ors_pulse_flipangle * DEG_RAD);

	for (int c = 0; c < N_protons; c++)
	{
		vec3& magn = protons[c].magnetization;

		float B_t = B_tot(protons[c]);
		float offset = gamma * B_t - gamma; // relative larmor frequency

		std::cout << offset << std::endl;

		protons[c].magnetization = (rot * vec4(protons[c].magnetization, 1.0f)).get_xyz();
	}
}

float Simulator::B_tot(Proton& p_proton)
{
	float sum = 0.0f;

	for (int i = 0; i < N_particles; i++)
		sum += B_dip(p_proton.position - particles[i].position);

	return B0 + sum;
}

float Simulator::B_dip(const vec3& p_r)
{
	float theta = vec3(0, 0, 1).angle(p_r);
	return B_eq * pow(particle_radius / p_r.length(), 3) * (3 * pow(cos(theta), 2) - 1.0f);
}