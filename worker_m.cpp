#include <thread>
#include <cassert>
#include <iostream>

#include "worker_m.hpp"

void Worker::write_data(size_t i, std::vector<std::vector<double> > d) {
//	std::cout << i << " write data" << std::endl;
	for (size_t m = 0; m < data.size(); ++m)
		data[m][i][i] = std::vector<double>(d[m]);
}

void Worker::read_data(size_t i) {
//	std::cout << i << " read data" << std::endl;
	for (size_t m = 0; m < data.size(); ++m)
		for (size_t j = 0; j < i; ++j) {
			data[m][i][j] = std::vector<double>(data[m][i - 1][j]);
		}
}

Worker::Worker(size_t id): BasicWorker(id) { }
void Worker::write(std::vector<std::vector<double>> d) { write_data(id, d); }
void Worker::read() { read_data(id); }

// for testing purposes
void Worker::run() {
	for (size_t i = 0; i < num_iterations; ++i) {
		if (id > 0)
			scheduler.wait_for_data_and_read(this);
		scheduler.send_data_to_worker(this, std::vector<std::vector<double> > {{10 * counter + (double)id}});
		if (id < num_workers - 1)
			scheduler.wait_for_worker_to_read(this);
		++counter;
	}
}
