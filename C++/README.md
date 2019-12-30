# MPCC
This is a C++ implementation of the MPCC controller. The implementation is NOT the version which was used in the paper [Optimization‐based autonomous racing of 1:43 scale RC cars](https://onlinelibrary.wiley.com/doi/abs/10.1002/oca.2123) but a new version with some additional features inspired the work we did in this paper [AMZ Driverless: The Full Autonomous Racing System](https://arxiv.org/abs/1905.05150)

## Difference to Matlab implementation

### Solver
This version only supports hpipm as a solver. However, the goal is to add an [acados](https://github.com/acados/acados) qp interface, which would allow to use the large list of solvers supported by acoados.

### Input rate cost and constraints
Instead of lifting the state and using the difference in inputs as a new input, this version uses the continuous time approach where the new inputs are the rate of change of the inputs. For more details look at the formulation in [AMZ Driverless: The Full Autonomous Racing System](https://arxiv.org/abs/1905.05150).

### Tire constraints
The C++ implementation adds the tire constraints used in [AMZ Driverless: The Full Autonomous Racing System](https://arxiv.org/abs/1905.05150). More precisely, I added a slip angle constraint (of the form: alpha_min <= alpha_f <= alpha_max) for the front wheel (since the 1:43 scale cars are rear wheel drive and have no brakes), and a tire friction ellipse constraint for the rear wheel (which has the form (e_long F_rx)^2 + F_ry^2 <= e_eps D_r). Note that if the car is all wheel drive or has breaks at the front wheel, also a tire ellipse constraint should be used for the front tire. Furthermore, it would be good practice to scale the tire friction forces, but for the cars used the forces are around 1, so this is not yet done.

### Obstacle Avoidance
There is no obstacle avoidance available yet in the C++ version 

## Installation 

To install all the dependencies run
```
./install.sh

```
this clones `blasfeo`, `hpipm`, `matplotlip-cpp`, `nlohmann/json`, and `eigen`, from their git repo, and safes them in a folder External. Additionally it installs `blasfeo` and `hpipm` in the same External folder, thus no admin rights are necessary.

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

There are still several things which should be added to the project. Most of them are marked with TODO in the code

