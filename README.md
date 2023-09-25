## README.md

### GeoRadarWithJulia

This repository contains a 2D Finite-Difference Time-Domain (FDTD) simulation of electromagnetic wave propagation with Perfectly Matched Layer (PML) absorbing boundary conditions. The code is implemented using Julia and is designed to model ground penetrating radar scenarios.

#### Background
The FDTD method is a numerical technique used to solve Maxwell's equations in a time-stepped manner. PMLs are designed to minimize reflections from computational domain boundaries, ensuring an accurate simulation of wave propagation in unbounded domains. The code presented in this repository is based on the method described in "Perfectly Matched Layer for the Absorption of Electromagnetic Waves" by Jean-Pierre Bérenger, 1994.

#### Repository Content

- 2D FDTD simulation with PML boundary conditions
- Source, field, and material structures for clear and organized data handling
- Time-stepping and update functions for the electromagnetic fields
- Utility functions for visualizing the simulation in real-time and recording results

#### Code Structure

1. Constants and Parameters: The code starts by defining various constants and parameters like the speed of light, impedance of free space, wavelength, etc.
2. Structures: Structures like `Field`, `Material`, and `Grid` are defined to keep the data organized.
3. Initialization Functions: These functions (`init_field`, `init_material`, and `init_grid_and_pml`) initialize the electromagnetic fields, material properties, grid, and PML parameters.
4. Field Update Functions: The `update_field!` function updates the electromagnetic fields using the FDTD method and considering the PML conditions.
5. Source Functions: The `gaussian_source` function returns a Gaussian pulse for use as the electromagnetic source in the simulation.
6. Visualization Function: The `plot_loop_field` function provides real-time visualization of the electric field.
7. Main Function: The `main` function initializes the necessary parameters, runs the simulation for the desired number of time steps, and visualizes the results.

#### How to Run

To execute the simulation:
1. Ensure you have Julia and the `PyPlot` library installed.
2. Download or clone the repository.
3. Run the provided Julia script.

#### Results

The simulation will visualize the propagation of the electromagnetic wave in real-time. At the end of the simulation, a recorded trace of the electric field at a specified location will be saved as "trace_rec.png".

#### References

Bérenger, Jean-Pierre. "A perfectly matched layer for the absorption of electromagnetic waves." Journal of computational physics 114.2 (1994): 185-200.

#### Acknowledgments

This code was inspired by various tutorials and academic references on FDTD simulations and PML implementations.

#### License

This project is licensed under the ETH License.

#### Contact

For questions, suggestions, or contributions, please open an issue or pull request on GitHub.
