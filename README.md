# MPCC
Simulation environment of the Model Predictive Contouring Controller (MPCC) for Autonomous Racing developed by the Automatic Control Lab (IfA) at ETH Zurich

## Fromulation

The MPCC is a model predictive path following controller which does follow a predefined reference path X^ref and Y^ref. This is achieved by augmenting the system with an integrator stated theta which approximates the progress along the reference path. The theta state is coupled to the real dynamics using the lag error e^l which is penalized in the cost. Additionally, the contouring error (lateral error to the reference path) is also penalized in the cost function. Finally, the progress along the reference path is maximized to achieve that the car does follow the path as fast as possible and the rate of the inputs is penalized. To guarantee that the car stays within the track, track constraints are imposed as well as bounds on all states, inputs and the rate of inputs. The resulting optimization problem is shown in the following equation:
<img src="https://github.com/alexliniger/MPCC/blob/master/Images/NMPCC_problem.jpg" width="700" />

The vehicle dynamics considered is a bicycle model with nonlinear magic formula tire models:
<img src="https://github.com/alexliniger/MPCC/blob/master/Images/Model.jpg" width="700" />

with the tire model and drive train model given as follows:
<img src="https://github.com/alexliniger/MPCC/blob/master/Images/TireModel.jpg" width="700" />

Finally, the state and inputs of the problem are given as follows:
<img src="https://github.com/alexliniger/MPCC/blob/master/Images/state-input.jpg" width="700" />

Where (X,Y) is the global position phi the heading of the car, v_x and v_y the longitudinal respectively the lateral velocity and omega the yaw rate. Theta is the augmented state which approximates the progress. The inputs are the duty cycle d to the drive train, the steering angle delta and the velocity along the reference path v_theta (which is an approximation of the velocities projected onto the path)

To achieve obstacle avoidance of other cars, before solving the MPCC problem grid search 1-D dynamic programming approach finds the best way to avoid the obstacles. This path is then converted into a corridor which modifies the track constraint. 

Finally, to solve the MPCC problem, we approximate the nonlinear problem as a time-varying quadratic program, by linearizing and discretizing the dynamics, approximating the lag and contouring errors using a first order tyler approximation and approximating the track constraints with two linear half space aligned with the boundaries. The resulting quadratic program is then solved using [hpipm](https://github.com/giaf/hpipm), however we also offer a Yalmip and CVX interface.

For more details see our paper [Optimization‐based autonomous racing of 1:43 scale RC cars](https://onlinelibrary.wiley.com/doi/abs/10.1002/oca.2123) or the [Arxiv](https://arxiv.org/abs/1711.07300) version.

## How to run

### Before running code
1) Install [hpipm](https://github.com/giaf/hpipm) including the python interface
2) alternativly install Yalmip or CVX
### Run code
0) in simulation.m change to the optimization framework you use (hpipm, Yalmip, CVX)
1) run simulation.m
2) play with the tunning in getMPC_vars.m
3) change the car model between FullSize and ORCA
4) change the track layout between the ORCA and the RCP track

### Notes
1) If you use the RCP track the obstacle postions need to be changed by hand

## Options

Beside the 1:43 RC cars, we also implemented a full sized car model as well as two possible tracks.

## Example
<img src="https://github.com/alexliniger/MPCC/blob/master/Images/MPC_sim.gif" width="700" />
