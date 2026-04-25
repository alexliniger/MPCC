# MPCC
This is a C++ implementation of the MPCC controller. The implementation is NOT the version which was used in the paper [Optimization‐based autonomous racing of 1:43 scale RC cars](https://onlinelibrary.wiley.com/doi/abs/10.1002/oca.2123) but a new version with some additional features inspired the work we did in this paper [AMZ Driverless: The Full Autonomous Racing System](https://arxiv.org/abs/1905.05150)

## Use With Full-Sized Car Model
The code in the master branch does not work well with the parameters of a full-sized car. However, I added a new branch (full-size) that uses a different SQP formulation which works for model parameters of a full-sized car. The simulation in this branch is performed for a formula student style car on a full-sized race track, following a pre-computed ideal line.

## Difference to Matlab implementation

### Solver
This version only supports hpipm as a solver. However, the goal is to add an [acados](https://github.com/acados/acados) qp interface, which would allow us to use the large list of solvers supported by acoados.

### Input rate cost and constraints - Dynamics
Instead of lifting the state and using the difference in inputs as new inputs, this version uses the continuous time approach where the new inputs are the rate of change of the inputs, similar to [AMZ Driverless: The Full Autonomous Racing System](https://arxiv.org/abs/1905.05150).

In detail we use the following dynamics,
<img src="https://github.com/alexliniger/MPCC/blob/master/Images/model_cpp.jpg" width="700" />

this also includes some changes in the notation, to match better the literature. Mainly the yaw rate is not `r` and the progress `s`. Thus, we have the following states and inputs,

<img src="https://github.com/alexliniger/MPCC/blob/master/Images/state_input_cpp.jpg" width="700" />

We also split up the force in x-direction into two components, the force at the wheel `F_r,x` and the friction force `F_fric`, which are defined as follows,
<img src="https://github.com/alexliniger/MPCC/blob/master/Images/forces_cpp.jpg" width="700" />

### Tire constraints
The C++ implementation adds the tire constraints used in [AMZ Driverless: The Full Autonomous Racing System](https://arxiv.org/abs/1905.05150). More precisely, I added a slip angle constraint for the front wheel (since the 1:43 scale cars are rear wheel drive and have no brakes), and a tire friction ellipse constraint for the rear wheel. Thus, the MPC problem the following three constraints, on top of state and input bounds,
<img src="https://github.com/alexliniger/MPCC/blob/master/Images/constraints_cpp.jpg" width="700" />

Note that if the car is all wheel drive or has brakes at the front wheel, also a tire ellipse constraint should be used for the front tire. 

### Beta Cost
We added an additional regularization cost, which penalizes high sideslip angles. This second regularization cost augments the small cost on the yaw rate. These regularization costs force the car to behave more reasonably and help the NLP to converge better.
<img src="https://github.com/alexliniger/MPCC/blob/master/Images/cost_cpp.jpg" width="700" />

### Obstacle Avoidance
There is no obstacle avoidance available yet in the C++ version 

### Options
Currently, only one track and car model is implemented. However, adapting the parameters only requires changing the json files in the Params folder.

## Installation 

To install all the dependencies run
```
./install.sh
```
this clones `blasfeo`, `hpipm`, `matplotlip-cpp`, `nlohmann/json`, and `eigen`, from their git repo, and safes them in a folder External. Additionally, it installs `blasfeo` and `hpipm` in the same External folder, thus no admin rights are necessary.

Note that `matplotlib-cpp` does also require `Python-2.7` and `matplotlib`, for more details see (https://github.com/lava/matplotlib-cpp).

Once all dependencies are installed `cmake` can be used to build the project
```
cmake CMakeLists.txt
make
```
To run the code simply execute the `MPCC`
```
./MPCC
```

### TODO

There are still several things that should be added to the project. Most of them are marked with TODO in the code.
