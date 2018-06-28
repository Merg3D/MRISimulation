#include "Simulator.h"

#include <fstream>

#include "Mat4.h"


Simulator::Simulator()
{
	// set all members
	max_iterations = 2e2;
	
	D = 3.0e-9;								// m^2 s^-1

	double T1 = 3.0;						// s
	double T2 = 80.0e-3;					// s

	R1 = 1.0 / T1;							// Hz
	R2 = 1.0 / T2;							// Hz

	Ms = 100.0;								// Am^2 / kg [Fe]
	M0 = 1.0;

	double mu_0 = 4e-7 * PI;
	B_eq = 0.72 * 5185 * mu_0 / 3.0 * Ms;	// Tesla (0.156)

	particle_radius = 6.08e-9;				// m8

	gamma = 267.513e6;						// rad s^−1 T^−1)

	normal_dt = pow(particle_radius, 2) / (6 * D);	// s
	dt = 5.55e-5;								// s

	exc_pulse_flipangle = 90.0;				// degrees
	ors_pulse_flipangle = 90.0;				// degrees

	volume = vec3d(1.5e-5);					// m

	B0 = 7.0;								// Tesla
	B1 = 0.0;								// Tesla

	ors_frequency = 0.0;
	ors_bandwidth = 100.0;
	gradient_duration = 20e-3;
	exc_pulse_time = 9999 * dt;

	N_threads = 4;

	if (gradient_duration > exc_pulse_time)
		gradient_duration = exc_pulse_time;

	int c = 0;
	float step = 0.1f;

	// set up experiments
	for (double frequency = 600; frequency <= 600; frequency += 100.0)
	{
		for (double bandwidth = 600.0; bandwidth <= 1100.0; bandwidth += 100.0)
		{
			for (double Cc = step; Cc <= 2.0f; Cc += step)
			{
				int max_particles = 500;

				double scale = max_particles / 6.38e12;
				double V = 2.5 * scale / Cc;			// ml

				Experiment exp;
				exp.ors_frequency = frequency;
				exp.ors_bandwidth = bandwidth;
				exp.N_particles = max_particles;
				exp.volume = vec3d(std::pow(V, 1.0 / 3.0) * 1e-2); // m^3
				exp.Cc = Cc;
				exp.id = c;
				exp.dt = dt;
				experiments.push_back(exp);

				c++;
			}
		}
	}

	std::cout << "Starting " << experiments.size() << " experiments" << std::endl;
}

Simulator::~Simulator()
{
	protons.clear();
	particles.clear();
	
	for (int c = 0; c < threads.size(); c++)
		threads[c].join();

	threads.clear();
}

void Simulator::start()
{
	for (int c = 0; c < experiments.size(); c++)
		start_experiment(experiments[c]);
}

void Simulator::start_experiment(Experiment& p_exp)
{
	current_exp = p_exp.id;

	std::cout << std::endl;
	std::cout << std::endl;
	std::cout << "Starting experiment " << current_exp << std::endl;

	ors_frequency = p_exp.ors_frequency;
	ors_bandwidth = p_exp.ors_bandwidth;
	N_particles = p_exp.N_particles;
	N_protons = p_exp.N_protons;
	dt = p_exp.dt;
	max_iterations = p_exp.max_iterations;
	volume = p_exp.volume;

	for (int a = 0; a < p_exp.averages; a++)
	{
		init();

		apply_ors_pulse();

		/*while (true)
		{
			if (iteration >= max_iterations)
				break;

			get_signal();
			iterate();

			t += dt;
			iteration++;
		}*/

		apply_exc_pulse();
		p_exp.results.push_back(get_signal());

		save();

		std::cout << std::endl;
	}

	save();
}

void Simulator::init()
{
	protons.resize(N_protons);
	particles.resize(N_particles);
	offsets.clear();
	signals.clear();

	for (int c = 0; c < N_protons; c++)
	{
		protons[c].position = vec3d(random.get_random(), random.get_random(), random.get_random()) * volume;
		protons[c].magnetization = vec3d(0, 0, M0);
	}

	for (int c = 0; c < N_particles; c++)
	{
		particles[c].position = vec3d(random.get_random(), random.get_random(), random.get_random()) * volume;
	}

	t = 0.0;
	iteration = 0;
	has_applied_exc_pulse = false;
	applying_gradient = true;
	running = true;
}

void Simulator::iterate()
{
	step_size = sqrt(6 * D * dt);

	for (int t = 0; t < N_threads; t++)
	{
		auto function = &Simulator::update_proton_positions;
		int start = t * N_protons / N_threads;
		int end = (t + 1) * N_protons / N_threads;

		threads.push_back(std::thread(function, this, start, end));
	}

	for (int t = 0; t < N_threads; t++)
		threads[t].join();

	threads.clear();

	for (int t = 0; t < N_threads; t++)
	{
		auto function = &Simulator::update_proton_magnetizations;
		int start = t * N_protons / N_threads;
		int end = (t + 1) * N_protons / N_threads;

		threads.push_back(std::thread(function, this, start, end));
	}

	for (int t = 0; t < N_threads; t++)
		threads[t].join();

	threads.clear();
}

double Simulator::get_z_magnetization()
{
	double sum = 0.0;

	for (int c = 0; c < N_protons; c++)
		sum += protons[c].magnetization.z;

	return sum / (double) N_protons;
}

double Simulator::get_signal()
{
	calculate_global_magnetization();
	vec3d M = global_magnetization;
	double signal = sqrt(M.x * M.x + M.y * M.y);

	signals.push_back(signal);
	std::cout << "Iteration: " << iteration << ", time: " << t << ", signal = " << signal << ", z = " << get_z_magnetization() << "\r";

	return signal;
}

void Simulator::calculate_global_magnetization()
{
	vec3d sum = 0.0;

	for (int c = 0; c < N_protons; c++)
		sum += protons[c].magnetization;

	global_magnetization = sum / (double) N_protons;
}

void Simulator::assert_in_space(Proton& p_proton)
{
	if (p_proton.position.x > volume.x)
		p_proton.position.x -= volume.x;
	else if (p_proton.position.x < 0.0)
		p_proton.position.x += volume.x;

	if (p_proton.position.y > volume.y)
		p_proton.position.y -= volume.y;
	else if (p_proton.position.y < 0.0)
		p_proton.position.y += volume.y;

	if (p_proton.position.z > volume.z)
		p_proton.position.z -= volume.z;
	else if (p_proton.position.z < 0.0)
		p_proton.position.z += volume.z;
}

void Simulator::update_proton_positions(int start, int end)
{
	for (int c = start; c < end; c++)
	{
		Proton& proton = protons[c];

		update_proton_position(proton);
		assert_in_space(proton);
	}
}

void Simulator::update_proton_magnetizations(int start, int end)
{
	for (int c = start; c < end; c++)
		update_proton_magnetization(protons[c]);
}

void Simulator::update_proton_position(Proton& p_proton)
{
	p_proton.position += random.get_random_vec3() * step_size;
}

void Simulator::update_proton_magnetization(Proton& p_proton)
{
	vec3d& magn = p_proton.magnetization;

	double Bt = B_tot(p_proton);
	double omega = gamma * Bt;

	//if (applying_gradient)
	//	omega += gamma * B1;

	double offset = omega - gamma * B0; // relative larmor frequency
	double theta = omega * dt;

	double E1 = exp(-R1 * dt);
	double E2 = exp(-R2 * dt);

	magn.x = E2 * (magn.x * cos(theta) + magn.y * sin(theta));
	magn.y = E2 * (magn.x * sin(theta) - magn.y * cos(theta));
	magn.z = E1 * magn.z + (1.0 - E1) * M0;
}

void Simulator::apply_exc_pulse()
{
	std::cout << std::endl;
	std::cout << "Applying exc pulse" << "\r";

	mat4 rot;
	rot.rotate_x(exc_pulse_flipangle * DEG_RAD);

	for (int c = 0; c < N_protons; c++)
		protons[c].magnetization = (rot * vec4d(protons[c].magnetization, 1.0)).get_xyz();
	
	std::cout << "Finished exc pulse" << std::endl;
}

void Simulator::apply_ors_pulse()
{
	mat4 rot;
	rot.rotate_x(ors_pulse_flipangle * DEG_RAD);

	std::cout << std::endl;
	std::cout << "Applying ORS pulse" << "\r";

	int count = 0;

	for (int c = 0; c < N_protons; c++)
	{
		vec3d& magn = protons[c].magnetization;

		double B_t = B_tot(protons[c]);
		double offset = gamma * B_t - gamma * B0; // relative larmor frequency

		if (offset > ors_frequency - ors_bandwidth / 2.0 && offset < ors_frequency + ors_bandwidth / 2.0)
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

void Simulator::save()
{
	std::ofstream file_off;
	std::ofstream file_decay;
	std::ofstream file_sig;
	std::string exp = "0";

	bool save_offsets = false;
	bool save_decay = true;

	if (save_offsets)
	{
		file_off.open("H:\\MyDocs\\BachelorProject\\Simulations\\data\\offsets.txt");

		for (int c = 0; c < offsets.size(); c++)
		{
			file_off << offsets[c] << std::endl;
		}

		file_off.close();
	}

	if (save_decay)
	{
		file_decay.open("H:\\MyDocs\\BachelorProject\\Simulations\\data\\decay.txt");

		for (int c = 0; c < signals.size(); c++)
		{
			file_decay << signals[c] << std::endl;
		}

		file_decay.close();
	}

	file_sig.open("H:\\MyDocs\\BachelorProject\\Simulations\\data\\data.txt");
	file_sig << "dt, volume, N_protons, N_particles, ors_bandwidth, ors_frequency, Cc, values" << "\n";
	
	for (int c = 0; c < experiments.size(); c++)
	{
		std::string result = "";

		for (int a = 0; a < experiments[c].results.size(); a++)
		{
			result += std::to_string(experiments[c].results[a]) + (a != experiments[c].results.size() -1 ? ", " : "");
		}

		int l = experiments[0].results.size() - experiments[c].results.size();

		if (experiments[c].results.size() > 0 && c != 0 && l != 0)
			result += ", ";

		for (int a = 0; a < l; a++)
		{
			result += std::string("1.000000") + (a != l - 1 ? ", " : "");
		}

		file_sig <<
			experiments[c].dt << ", " <<
			experiments[c].volume.x << ", " <<
			experiments[c].N_protons << ", " <<
			experiments[c].N_particles << ", " <<
			experiments[c].ors_bandwidth << ", " <<
			experiments[c].ors_frequency << ", " <<
			experiments[c].Cc << ", " <<
			result << "\n";
	}
	
	file_sig.close();
}
