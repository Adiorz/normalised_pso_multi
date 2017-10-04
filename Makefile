main_m: additional.cpp  data_preprosessor.cpp  fft.cpp  main_m.cpp  pso_m.cpp  scheduler_m.cpp  worker_m.cpp
	g++ additional.cpp  data_preprosessor.cpp  fft.cpp  main_m.cpp  pso_m.cpp  scheduler_m.cpp  worker_m.cpp -std=c++14 -pedantic -lpthread -o main_m -ggdb

clean:
	rm -f main_m
