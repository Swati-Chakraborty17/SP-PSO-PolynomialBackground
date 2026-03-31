# PSO-Based Self-Potential Optimization Using a Dipping Thin Sheet Model with Modeled quadratic 
# background effect

This repository implements a **standard Particle Swarm Optimization (PSO) workflow**
with **explicitly tested control parameters** for the inversion of **self-potential (SP) data**
using a **dipping thin-sheet source model** combined with a **quadratic regional background**.

The code is designed for **research and reproducibility**, with multiple independent PSO
runs used to assess solution stability and parameter uncertainty.
NOTE: Since PSO is an optimizer and not a deterministic sampling method, solution stability is
evaluated using ensemble-inspired statistics derived from multiple independent optimization runs.

---

## Repository Contents

### Main Program
- **sp_pso_optimization.m**  
  Executes the complete inversion workflow, including data loading, multi-run PSO
  optimization, convergence analysis, ensemble statistics, and visualization.

### Core Functions
- **pso.m**  
  Particle Swarm Optimization algorithm with linearly decreasing inertia weight,
  boundary reflection handling, and full history tracking of particle positions
  and fitness values.

- **model_thin_sheet.m**  
  Forward model for the self-potential anomaly of a dipping thin sheet, including
  a quadratic regional background term of the form:
  
   V_{bg}(x) = c + b x + i x^2
 
- **objective_function.m**  
  Objective (misfit) function based on percentage RMS error, incorporating a soft
  boundary penalty to discourage boundary-hugging solutions while maintaining
  smooth optimization behavior.

### Data
- **Model_eta.dat**  
  Input data file containing:
  - Column 1: Profile distance (m)
  - Column 2: Observed self-potential anomaly (mV)

---

## How to Run

1. Place all `.m` files and `Model_eta.dat` in the same directory.
2. Open MATLAB and set the working directory to this folder.
3. Run the main script:   'sp_pso_optimization.m'
