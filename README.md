# swPHoToNs
	swPHoToNs(Sunway Parallel Heterogeneous & Threads oriented code for cosmological N-body simulation) is a Particle-Mesh(PM) and Fast Multipole Method(FMM) based code that can perform cosmological simulations with trillions of particles efficiently on the Sunway TaihuLight supercomputer. swPHoToNs is written in C and Fortran. We managed to run the original version and optimized version of swPHoToNs on the Sunway TaihuLight supercomputer. The source code for prior versions is publicly available. We have performed a series of cosmological simulations with different problem sizes to evaluate the performance and scalability of our code. The initial condition of the simulation is setup by particle displacement, according to the Zel’dovich approximation method. The largest weak scaling experiment contains up to 1.6 trillion particles and is scaled to 10,400,000 cores with a parallel efficiency of 80.9% and a sustained performance of 56.3 PFlops.

## Compile
	Unzip the 2DECOMP&FFT package and make the library.
	Use "make" command to compile the source code and link object files.

## Run
	Enter "run" directory, run suh64.sh.


If there are any issues, you have questions, do not hesitate to send an email to <liuzhao@mail.nsccwx.cn>.
