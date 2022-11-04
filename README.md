# Arbor: Ornstein-Uhlenbeck on LIF example

This is an example to play around with, featuring a simple LIF neuron with Ornstein-Uhlenbeck input current, implemented in Arbor.

The file `set_arbor_env` will need to be adapted to the local Arbor installation.

Execute the `build` shell script to build the custom catalogues.

Execute the `run` shell script to run the example in Python.

Tested with Arbor version 0.6.1-dev, state of commit [8af6bd2](https://github.com/arbor-sim/arbor/commit/8af6bd273678406bab7bbd0abd5c29b2d91ba6cd) (including SDE computation as of commit [5d141aa](https://github.com/arbor-sim/arbor/pull/1884/commits/5d141aae533aecb81c00cf8b26925d5dff09478e)).
