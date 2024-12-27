# openAbsorptionSolver

A software package to calculate open absorption of falling film absorbers using OpenFOAM v7 (org).

---

## Description

This software package simulates the absorption of a liquid desiccant solution film falling down a vertical plate. Adaptations for different geometries are possible. The liquid desiccant solution is assumed to be a binary mixture of water and a desiccant, such as lithium bromide. 

The absorption process is driven by the vapor pressure difference between the liquid desiccant solution and the air flow, resulting in mass transfer of water vapor. This process involves heat transfer, including latent heat. A schematic of the geometry is shown below:

![Geometry](https://ars.els-cdn.com/content/image/1-s2.0-S1359431123002120-gr1.jpg)

### Key Features
1. **Custom OpenFOAM solvers:**  
   - `heatMassTransferFoam`  
   - `heatMassTransferWireFoam`
2. **Custom boundary functions:**  
   - Couple manager and region coupling utilities.
3. **Custom field functions:**  
   - Absorbed mass flux, heat flux, and related utilities.
4. **Python scripts for case generation:**  
   - `wire_create_case.py`  
   - `stacked_create_case.py`
5. **Supporting scripts for automation and analysis:**  
   - Parallelization, mapping, and post-processing utilities.

### Directory Structure
- `./baseCase`: Base OpenFOAM cases for stacked and wire geometries. 
- `./helpingScripts`: Additional supporting scripts for case setup and simulation.  
- `./pythonScripts`: Python scripts for automated case generation.  
- `./solvers`: Custom OpenFOAM solvers and libraries.  
- `./workingFolder`: Workspace for generated cases.  
- `./README.md`: This file.  
- `./LICENSE`: License file.

More information on usage and examples can be found in the **Usage Advice** section. A more detailed explanation of the technical details of the solver and a comparison with experimental data can be found in the open access publication by Hochwallner et al. (2023), listed in the **References** section. A future publication by Hochwallner et al. (2025), currently under review, will provide a more detailed explanation of the current state of the solver. 

---

## Installation Instructions

### Prerequisites
- **OpenFOAM v7 (org)** with `wmake` for compilation.  
- **Python with Conda** for running automation scripts.

### Steps

#### 1. Compile Custom Boundary Functions
Navigate to `./solvers/libraries/customBoundaries` and run:
```bash
wclean
wmake
```
Output: `$(FOAM_USER_LIBBIN)/customBoundaryConditions`.

#### 2. Compile Custom Field Functions
Navigate to `./solvers/libraries/customFieldFunctionObjects` and run:
```bash
wclean
wmake
```
Output: `$(FOAM_USER_LIBBIN)/customFieldFunctionObjects`.

#### 3. Compile Custom Solvers
Navigate to each solver directory (`./solvers`) and run:
```bash
wclean
wmake
```
Output: `$(FOAM_USER_APPBIN)`.

#### 4. Set Up Python Environment
Navigate to `./pythonScripts` and execute:
```bash
conda env create -f conda-env.yml
conda activate openAbsorptionSolver-env
```

---

## Usage Advice

The software was intended as a tool to optimize the absorption of falling film absorbers with vertical plates. For this reason two different geometries are implemented, a stacked arrangement and a wire case. The software is designed to be flexible and adaptable to other geometries. 

### Overview of Geometries
1. **Stacked Arrangement:** Consecutive offset falling films to bring the humid air core into contact with the liquid desiccant film. 

2. **Wire Geometry:** Multiple wires in the air gap inducing a Kármán vortex street to enhance heat and mass transfer.  

(Pictures will be added as soon as the current publication is published.)

### Simulation Process

#### Stacked Arrangement
1. **Generation of initial case**  
   This case is a 2D case with a single cell in the z-direction and is used to calculate the steady air flow between the falling films. The case is generated using the Python script `stacked_create_init_case.py`.

2. **Optional: Decomposition of the initial case**  
   The initial case can be decomposed using the OpenFOAM utility `decomposePar`. This is automatically done by the Python script, if set.

3. **Simulation of the initial case**  
   The initial case is simulated to reach a steady state. The simulation is done using the OpenFOAM solver `simpleFoam`.

4. **Check of the initial case**  
   The initial case is checked for convergence, and the results are analyzed for physical plausibility.

5. **Transformation of the initial case**  
   The 2D mesh of the initial case is transformed to the height of the falling film absorber. The transformation is done using the OpenFOAM utility `transformPoints -scale "(1 1 X)"`, where X is the ratio of the height of the falling film absorber to the height of the initial case. For the given example case, `X` is 5000.

6. **Optional: Reconstruction of the initial case**  
   If the initial case was parallelized, the case has to be reconstructed using the OpenFOAM utility `reconstructPar`.

7. **Generation of the full case**  
   The case is generated using the Python script `stacked_create_case.py`. The flow solution of the solution film is mapped by the Python script, assuming a fully developed Nusselt profile for the full film.

8. **Mapping of air flow solution of the initial case on the air region of the full case**  
   The air flow solution of the initial case is mapped on the air region of the full case using the OpenFOAM utility `mapFields -consistent -targetRegion air -sourceTime X ../stacked_init`, where `X` is the final time of the initial case and `../stacked_init` is the directory of the initial case.

9. **Optional: Decomposition of the full case**  
   The full case can be decomposed to increase computational speed. The initial decomposition at the creation of the case by the Python script `stacked_create_case.py` is not valid anymore, as the air flow solution of the initial case is mapped on the air region of the full case. However, the initial decomposition is necessary to compute the necessary cell allocations for this decomposition. The decomposition is done using the OpenFOAM utilities `decomposePar -region air`, `decomposePar -region solution1`, and `decomposePar -region solution2`.

10. **Simulation of the full case**  
    The full case is simulated using the OpenFOAM solver `heatMassTransferFoam`.

11. **Optional: Reconstruction of the full case**  
    If the full case was parallelized, the case has to be reconstructed using the OpenFOAM utilities `reconstructPar -region air`, `reconstructPar -region solution1`, and `reconstructPar -region solution2`.

#### Wire Geometry
1. **Generation of initial case**  
   This case is a 2D case with a single cell in the z-direction and is used to calculate the transient air flow between the wires. The case is generated using the Python script `wire_create_init_case.py`.

2. **Optional: Decomposition of the initial case**  
   The initial case can be decomposed using the OpenFOAM utility `decomposePar`. This is automatically done by the Python script, if set.

3. **Simulation of the initial case**  
   The initial case is simulated using the OpenFOAM solver `pimpleFoam`. The simulation time has to be chosen long enough to reach a quasi-steady state. For the given example case, this first simulation time was set to 2 seconds. After reaching the quasi-steady state, the simulation flow results have to be written in a sufficiently fine manner to depict the periodic flow. For the given example case, the written time steps were set to 1e-5 seconds.

4. **Optional: Reconstruction of the initial case**  
   If the initial case was parallelized, the case has to be reconstructed using the OpenFOAM utility `reconstructPar`.

5. **Check of the initial case**  
   The initial case is checked for convergence, and the results are analyzed for physical plausibility. Furthermore, the written time steps have to be checked to ensure they sufficiently depict the periodic flow. One way to check this is to plot the pressure drop over time and verify periodicity. If the solution is periodic, the period is calculated, and only the time steps of one period are used for the next simulation. Later time steps are deleted. This is done using the Python script `periodic_solution.py` in the `./helpingScripts` directory.

6. **Transformation of the initial case**  
   The 2D mesh of the initial case is transformed to the height of the full wire case. The transformation is done using the OpenFOAM utility `transformPoints -scale "(1 1 X)"`, where X is the ratio of the height of the wire case to the height of the initial case. For the given example case, `X` is 5000.

7. **Generation of the full case**  
   The case is generated using the Python script `wire_create_case.py`. The flow solution of the solution film is mapped by the Python script, assuming a fully developed Nusselt profile for the full film.

8. **Mapping of air flow solution of the initial case on the air region of the full case**  
   The solution of the initial case is mapped on the air region of the full case using the OpenFOAM utility `mapFields -consistent -targetRegion air -sourceTime X ../wire_init`, where `X` is the time step of the initial case and `../wire_init` is the directory of the initial case. This has to be done for all time steps of the initial case. The mapping is done using the bash script `map_fields.sh` in the `./helpingScripts` directory.

9. **Prepare the full case for the `heatMassTransferWireFoam` solver**  
   The `periodicFlow` directory and `periodicFlowDict` file have to be created to use the `heatMassTransferWireFoam` solver. In this folder, the time steps of the periodic air flow are written, and the time steps are listed in the `periodicFlowDict` file. This is done using the bash script `prepare_case.sh` in the `./helpingScripts` directory.

10. **Optional: Presimulate a steady-state solution**  
    To speed up the simulation, it is beneficial to start with a steady-state solution of the absorption problem before starting the transient simulation. A steady-state solution can be calculated using the `heatMassTransferFoam` solver. The bash script `simulate_steady_state.sh` in the `./helpingScripts` directory helps prepare the case for this, by taking the first time directory of the `periodicFlow` folder, modifying the `controlDict` and `fvSchemes` for steady-state, and adding relaxation factors. The steady-state solution is then calculated using the `heatMassTransferFoam` solver. To post-process the steady-state solution, the bash script `post_process_steady_state.sh` in the `./helpingScripts` directory can be used. The script writes the steady-state solution as the initial condition for the following transient simulation and reverts the changes in the `controlDict`, `fvSchemes`, and relaxation factors for the transient simulation.

11. **Optional: Decomposition of the full case**  
    The full case can be decomposed to increase computational speed. The initial decomposition at the creation of the case by the Python script `wire_create_case.py` is not valid anymore, as the air flow solution of the initial case is mapped on the air region of the full case. However, the initial decomposition is necessary to compute the required cell allocations for this decomposition. The decomposition is done using the bash script `decompose_wire.sh` in the `./helpingScripts` directory, which decomposes all time steps of the `periodicFlow` folder for all regions and writes the decomposed time steps into the `periodicFlow` folder of each processor folder.

12. **Simulation of the full case**  
    The full case is simulated using the OpenFOAM solver `heatMassTransferWireFoam`.

13. **Optional: Reconstruction of the full case**  
    If the full case was parallelized, the case has to be reconstructed using the OpenFOAM utilities `reconstructPar -region air`, `reconstructPar -region solution1`, and `reconstructPar -region solution2`.

---

## References

- Hochwallner, F., Reichl, C., Emhofer, J. (2021). Frostprävention in Kühlhäusern. *47th Annual Meeting of the Deutscher Kaelte und Klimatechnischer Verein 2021*, Vol. 1, pp. 1014-1027.
- Hochwallner, F., Reichl, C., Emhofer, J. (2023). Reduced modeling of liquid desiccant falling film absorbers. *Applied Thermal Engineering, 225*, 120183. [https://doi.org/10.1016/j.applthermaleng.2023.120183](https://doi.org/10.1016/j.applthermaleng.2023.120183)
- Hochwallner, F., Stracke, C., Emhofer, J., Reichl, C., 2025. Optimizing the absorption of falling film absorbers with vertical plates utilizing an open-source newly developed CFD solver. *Applied Thermal Engineering* (submitted).

---

## Acknowledgments

This work was funded by the **Austrian Research Promotion Agency (FFG)** under Grant No. 874186. The solver has been extensively tested using the resources of the Vienna Scientific Cluster (VSC).

---

## License

The software is released under the **GNU General Public License v3 (GPL v3)**, ensuring it remains free and open-source for use, modification, and distribution. See `./LICENSE` for details.
