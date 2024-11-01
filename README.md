Dependencies: 
Seaborn-0.13.2
TinyGP-0.3.0
SciPy 1.13.1
NumPy 2.0.2
Gfortran 11.4.1

If you want the figure to immediately be displayed at the end of the calculation, you also need Okular.


The program is run as follows:

`python3 run_experiment.py <element> <l_distribution> <initial_parameter guess(es))`

where the l_distribution is either:
    0 for modified statistical
    2 for quadratic 
    4 for statistical

This code uses the Akylas cascade code [1] to calculate muonic X-ray transitions for a given element. Using experimental data for whole shell transitions, a least squares fit is used to find appropriate parameters for the initial l-distributions.

These parameters are then used to form a mean l-distribution, which can be used as the mean in a Gaussian Process(GP). Samples can then be drawn from this GP and inputted into the cascade code. As a result, an intensity distribution can be created, showing the relationship between l-distribution and final intensity.


[1] V. R. Akylas and P. Vogel, Muonic atom cascade program, Computer Physics Communications 15, 291 (1978).

