#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <vector>
#include <algorithm>

#include "data_preprosessor.hpp"
#include "scheduler_m.hpp"
#include "worker_m.hpp"
#include "pso_m.hpp"
#include "fft.hpp"
#include "additional.hpp"

#define num_dims 3
#define numofparticles 50

// Declarations of static members
Scheduler Worker::scheduler = Scheduler(0);
size_t Worker::num_workers;
size_t Worker::num_iterations;
std::vector<std::vector<std::vector<std::vector<double> > > > Worker::data;

int main(int argc, char *argv[]) {
	// get input arguments
	std::string file_name(argv[1]);
	std::cout << "Filename: " << file_name << std::endl;
	std::string log_file_name(argv[2]);
	std::cout << "Log filename: " << log_file_name << std::endl;
	size_t channel = atoi(argv[3]);
	std::cout << "Channel: " << channel << std::endl;

	size_t numofsamples_2 = atoi(argv[4]);
	std::cout << "numofsamples_2: " << numofsamples_2 << std::endl;

	int W = atoi(argv[5]);
	std::cout << "Number of workers (swarms): " << W << std::endl; //modes number
	int I = atoi(argv[6]);
	std::cout << "Number of iterations: " << I << std::endl; //number of iterations for each worker (swarm)

	size_t numofmodes_per_swarm = 1;

	// pre-process data
	// 0 freqs, 1 hammer, 2-16 15 channels
	//size_t numofsamples_2 = 501;
	size_t numofsamples = (numofsamples_2 - 1) * 2;
	DataStream freq_stream(file_name, 0, numofsamples_2, 0);
	DataStream channel_stream(file_name, channel, numofsamples_2, 0);

	// store captured frequencies and corresponding amplitudes in double vectors
	std::vector<double> freq(freq_stream.get());
	channel_stream.fill_first_elements(10, 0);    // zero first 10 values
	std::vector<double> amp(channel_stream.get());
	std::vector<double> amp_gauss(gaussian_filter(amp, 1));

	double fs = numofsamples; //5000 is real for the data
	std::vector<double> time(numofsamples);

	get_times(time, numofsamples, 1 / fs);

	// init min/max values
	std::vector<std::vector<double> > xmin(numofmodes_per_swarm);
	std::vector<std::vector<double> > xmax(numofmodes_per_swarm);
	for (size_t m = 0; m < numofmodes_per_swarm; ++m)
		init_minmax(xmin[m], xmax[m], num_dims, amp, fs);

	Worker::scheduler = Scheduler(W);
	Worker::num_workers = W;
	Worker::num_iterations = I;
	Worker::data = std::vector<std::vector<std::vector<std::vector<double> > > >(numofmodes_per_swarm);

	for (size_t m = 0; m < numofmodes_per_swarm; ++m) {
		Worker::data[m] = std::vector<std::vector<std::vector<double> > >(W);
		for (size_t i = 0; i < W; ++i) {
			Worker::data[m][i] = std::vector<std::vector<double> >(i + 1);
			for (size_t j = 0; j < i + 1; ++j) {
				Worker::data[m][i][j] = std::vector<double>(num_dims);
			}
		}

//		for (size_t i = 0; i < W; ++i)
//			std::cout << Worker::data[m][i].size() << std::endl;
	}

    std::vector<PSO*> workers = std::vector<PSO*>();
	std::vector<std::thread> t_workers(W);

	for (size_t i = 0; i < W; ++i) {
//		std::cout << "pushing " << i << std::endl;
		workers.push_back(
				new PSO(numofparticles, numofmodes_per_swarm, num_dims, xmin, xmax, &time, &amp,
						numofsamples, i));
	}

	for (size_t i = 0; i < W; ++i) {
//		std::cout << "going to run " << i << std::endl;
		t_workers[i] = std::thread(&PSO::run, std::ref(*workers[i]));
	}

	for (size_t i = 0; i < W; ++i)
		t_workers[i].join();

	std::cout << "Found parameters" << std::endl;
	std::cout << "----" << std::endl;
	for (size_t m = 0; m < numofmodes_per_swarm; ++m)
		for (size_t i = 0; i < workers.size(); ++i) {
			std::cout << "m[" << m << "]: " << workers[i]->get_id() << " ";
			for (std::vector<double> x: workers[i]->getgbest()) {
				for (double x_m: x)
					std::cout << x_m << " ";
			}
			std::cout << std::endl;
		}
	std::cout << "----" << std::endl;

	// Prepare data log for visualisation
//	for (size_t m = 0; m < numofmodes_per_swarm; ++m)
//		prepare_log_file_for_visualisation(log_file_name, W, Worker::data[m][W - 1], time, amp, amp_gauss, numofsamples);



	// store founds for second round
	std::vector<std::vector<double> > founds(Worker::data[0][W - 1]);

	// we don't need 'old' workers anymore
	for (auto w: workers)
		delete w;

	numofmodes_per_swarm = W;
	std::cout << W << std::endl;

	xmin = std::vector<std::vector<double> > (numofmodes_per_swarm);
	xmax = std::vector<std::vector<double> > (numofmodes_per_swarm);

	//find max
	std::vector<std::vector<double> > data_copy(numofmodes_per_swarm);

	size_t f;
	size_t f_l;
	size_t f_r;
	std::vector<size_t> idxL;
	std::vector<size_t> idxR;
	for (size_t m = 0; m < numofmodes_per_swarm; ++m) {
		// determine frequency ranges for each of modes
		f = founds[m][1];
		findMinimas(amp_gauss, 0, f, idxL);
		findMinimas(amp_gauss, f, numofsamples_2-1, idxR);
		f_l = idxL[idxL.size() - 1];
		f_r = idxR[0];

		std::cout << "helper " << m << ": <" << freq[f_l] << ", " << f << ", " << freq[f_r] << ">" << std::endl;

		size_t f_range = f_r - f_l;
		std::vector<double>::const_iterator first = amp.begin() + f_l;
		std::vector<double>::const_iterator last = amp.begin() + f_r;

		data_copy[m] = std::vector<double>(first, last);
		std::sort(data_copy[m].begin(), data_copy[m].end(), &abs_compare);
		double max_abs = 2.0 * data_copy[m][0];
		std::cout << "max: " << max_abs << std::endl;

		xmin[m] = std::vector<double> ({0.0, freq[f_l], 0.0001});
		xmax[m] = std::vector<double> ({0.23, freq[f_r], 0.2});
	}

	Worker::scheduler = Scheduler(0);
	Worker::num_workers = 1;
	Worker::num_iterations = I;
	Worker::data = std::vector<std::vector<std::vector<std::vector<double> > > >(numofmodes_per_swarm);
	PSO pso(numofparticles, numofmodes_per_swarm, num_dims, xmin, xmax, &time, &amp,
							numofsamples, 0, true);
	pso.run();

	std::cout << "Found parameters" << std::endl;
	std::cout << "----" << std::endl;
	for (std::vector<double> x: pso.getgbest()) {
		for (double x_m: x)
			std::cout << doubleToText(x_m) << " ";
	}
	std::cout << std::endl << "----" << std::endl;

	prepare_log_file_for_visualisation(log_file_name, W, pso.getgbest(), time, amp, amp_gauss, numofsamples);

	return 0;
}
