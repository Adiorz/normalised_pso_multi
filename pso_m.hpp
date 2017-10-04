#pragma once

#include <vector>
#include <iostream>

#include "worker_m.hpp"

class PSO: public Worker {
private:
	size_t id;
	double c1, c2;
	double w;

public:
	bool initialized = false;

	size_t numofparticles;
	size_t numofmodes;
	size_t numofdims;
	size_t numofsamples;
	size_t numofsamples_2;

	bool helper;

	std::vector<std::vector<std::vector<double> > > V;
	std::vector<std::vector<std::vector<double> > > X;
	std::vector<std::vector<double> > Xmax;
	std::vector<std::vector<double> > Xmin;
	std::vector<double> Vmax;
	std::vector<double> Vmin;

	std::vector<std::vector<std::vector<double> > > pbests;
	std::vector<double> pbestfits;
	std::vector<double> fitnesses;
	std::vector<double> gbestfit;
	std::vector<size_t> sorted_indices;

	std::vector<std::vector<double> > gbest;
	std::vector<std::vector<double> > bests;
	//TODO: continue below

	std::vector<double> *A;
	std::vector<double> *time;
	std::vector<double> A_gauss;
	std::vector<double> P;

	std::vector<std::vector<size_t>> to_skip;

	double init_param(size_t m, size_t j);

	void init_max_velocities();

	void initpopulation();

	void addparticle();

	void reposition_particle(size_t m, size_t i);

//	double fitnessfunc_singleparticle(std::vector<double> &p, size_t m);
	double fitnessfunc_singleparticle(std::vector<std::vector<double>> p);

	static double calc_response(
			std::vector<std::vector<double> > results, double t,
			size_t id = 0);

	void fitnessfunc(bool first = false);

	void calcgbest(bool first = false);

	void update();

	double update_velocity(double w, double X, double V, double V_min,
			double V_max, double gbest, double pbests, double c1, double c2);

	double update_position(double X, double V, double X_min, double X_max);

public:
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution;

	PSO(size_t numofparticles, size_t numofmodes, size_t numofdims,
			std::vector<std::vector<double> > Xmin,
			std::vector<std::vector<double> > Xmax,
			std::vector<double> *time,
			std::vector<double> *A, size_t numofsamples,
			size_t id, bool helper = false, double c1 = 2, double c2 = 2);

	void run();

	std::vector<std::vector<double> > getgbest() const;

	std::vector<double> getgbestfit();
};
