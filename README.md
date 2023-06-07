Quasar Time Dilation
====================

(c) 2022, 2023 Geraint F Lewis, B. J. Brewer

LICENSE: GNU General Public Licence v3

Dependencies
============

`NumPy`, `matplotlib`, `scipy`, `astropy`,
`DNest4` (compiled with C extensions on)

Running
=======

To run the MCMC, execute `python main.py`. This will load the data and launch
Diffusive Nested Sampling. Every 100 iterations, the posterior samples are
saved to `posterior_sample.txt`, and the marginal likelihood is printed to the
screen. The main parameter of interest, n (which describes the dependence of
the quasar timescale on the redshift), is the first column of this text file.
Example output is given below.

```
Data read and file closed
# Seeding random number generators. First seed = 1234.
# Generating 5 particles from the prior...done.
# Creating level 1 with log likelihood = -4316.61360484.
# Creating level 2 with log likelihood = -3205.8457954.
Sample 100.
log(Z) = -1318.8044390992486
Information = 5.568074435862172 nats.
Effective sample size = 1.000000021225817
```

log(Z) is the marginal likelihood. For more information about Diffusive Nested
Sampling, see [this paper](https://www.jstatsoft.org/index.php/jss/article/view/v086i07/v86i07.pdf),
but be aware that the version shipped with the paper contains a bug. It is
recommended that you install from the `master` branch of the
[DNest4 git repository](https://github.com/eggplantbren/DNest4).

Changing the model
==================

By default, the model sampled is the favoured one, with the parameter n fixed
to 1. To modify the model, follow the instructions in the comment in
the `prior_transform()` function in `main.py`:

```
def prior_transform(u):
    """Transforms the uniform random variable `u ~ Unif[0., 1.)`
    to the parameter of interest `x ~ Unif[1, 5)`."""

    x = 1. + 4*u 

    # This is n set to 1. Change the first term here to set it to another
    # value, or comment out this line to let n be free.
    x[0] = 1 + 1e-6*(u[0]-0.5)

    return x
```


