# MPCC
This is a Python implementation of the MPCC controller. As of 27 January 2021, the implentation is under review. Please note that the controller does not work as a correct manner as Matlab/C++ implementation.

## Instllation
This implementation depends on the following packages. Please see `requirements.txt` for more details.

- matplotlib
- numpy 
- pandas
- scipy
- qpsolvers
- osqp
- numba
- sympy

## How to run
With `--render` option, you can watch the animation

```
python nmpc_sim.py --render
```
