# Klein-Gordon Field Simulator
This project presents a first-principles simulation of particle dynamics governed by the 1D
Klein-Gordon equation, a foundational model in quantum field theory. The theoretical framework involves discretizing the continuous field onto a spatial lattice and mapping the resulting
Hamiltonian to an equivalent XY model in a transverse field. To execute the simulation, a
high-performance statevector simulator was developed from scratch in C++ using the Eigen
library for efficient linear algebra. The simulator’s output is piped in real-time to a separate Python script, which leverages Matplotlib and SciPy to generate a smooth, animated
visualization of the particle’s probability wavepacket. The simulation successfully models
the time evolution of a single-particle excitation, clearly demonstrating the characteristic
propagation and delocalization of its wavepacket across the lattice, in agreement with theoretical predictions. This work serves as a comprehensive demonstration of the pipeline from
abstract field theory to concrete computational physics, providing a robust foundation for
future explorations into more complex interacting theories.

