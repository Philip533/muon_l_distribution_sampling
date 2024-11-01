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
