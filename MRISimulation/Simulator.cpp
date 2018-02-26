#include "Simulator.h"



Simulator::Simulator()
{
	D = 3.0e-9f; // m^2 s^-1


}


Simulator::~Simulator()
{
}

void Simulator::iterate()
{
	step_size = std::sqrt(6 * D * dt);

	for (int c = 0; c < particles.size(); c++)
	{
		Particle* p = particles[c];

		update_particle(p);
		assert_in_space(p);

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
		Proton* p = reinterpret_cast<Proton*>(particles[c]);
		sum += p->magnetization;
	}

	global_magnetization = sum / (float)particles.size();
}

void Simulator::assert_in_space(Particle* p_particle)
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

void Simulator::update_particle(Particle* p_particle)
{

}

void Simulator::update_particle_position(Particle* p_particle)
{

}

void Simulator::update_particle_magnetization(Particle* p_particle)
{
	Proton* p = reinterpret_cast<Proton*>(p_particle);
	vec3& magn = p->magnetization;

	float theta;
	float E1 = std::exp(-R2 * dt);
	float E2 = std::exp(-R1 * dt);

	magn.x += E1 * (magn.x * std::cos(theta) + magn.y * std::sin(theta));
	magn.y += E1 * (magn.x * std::sin(theta) - magn.y * std::cos(theta));
	magn.z += E2 * magn.z + (1.0f - std::exp(-dt / T1));
}
