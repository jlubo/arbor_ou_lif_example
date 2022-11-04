: Input current modeled by an Ornstein-Uhlenbeck process (given number and firing rate of neurons in a putative population)

NEURON {
	POINT_PROCESS ornstein_uhlenbeck_pop
	RANGE N, f, w_out, tau
	NONSPECIFIC_CURRENT I
}

UNITS {
	(ms) = (milliseconds)
	(nA) = (nanoampere)
}

PARAMETER {
	N = 10 : number of neurons in the putative population
	f = 20 (Hz) : mean firing rate of the neurons in the putative population
	w_out = 1 (nC) : synaptic weight of synapses outgoing from the putative poulation
	tau = 5.0 (ms) : time constant (basically synaptic time constant)
}

STATE {
	I_ou (nA) : instantaneous state of the OU process
	active : indicates if the process is currently enabled 
	mean (nA) : mean of the stochastic process
	stdev : standard deviation of the stochastic process (in units of nA*s^(1/2))
}

INITIAL {
	I_ou = 0
	active = -1
	mean = N * f * w_out : mean of the stochastic process
	stdev = (N * f)^(1/2) * w_out : standard deviation of the stochastic process (in units of nA*s^(1/2))
}

BREAKPOINT {
	SOLVE state METHOD stochastic
	
	I = -I_ou
}

WHITE_NOISE {
    zeta
}

DERIVATIVE state {
	: Stochastic process
	I_ou' = heaviside(active) * (-I_ou + mean + stdev * 1000^(1/2) * zeta) / tau
}

NET_RECEIVE(weight) {
	
	if (weight >= 0) { : indicates that stimulation begins
		I_ou = mean : initialize the process at the mean value
		active = 1
	}
	else { : indicates that stimulation ends
		I_ou = 0 : switch off the process
		active = -1
	}
}


FUNCTION heaviside(x) { : the Heaviside Theta function
	if (x >= 0) {
		heaviside = 1
	}
	else {
		heaviside = 0
	}
}
