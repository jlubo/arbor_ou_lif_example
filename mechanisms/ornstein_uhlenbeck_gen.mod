: Input current modeled by an Ornstein-Uhlenbeck process (given mean and standard deviation)

NEURON {
	POINT_PROCESS ornstein_uhlenbeck_gen
	RANGE mean, stdev, tau
	NONSPECIFIC_CURRENT I
}

UNITS {
	(ms) = (milliseconds)
	(nA) = (nanoampere)
}

PARAMETER {
	mean = 1 (nA) : mean of the stochastic process
	stdev = 1 : standard deviation of the stochastic process (in units of nA*s^(1/2))
	tau = 5.0 (ms) : time constant (basically synaptic time constant)
}

STATE {
	I_ou (nA) : instantaneous state of the OU process
	active : indicates if the process is currently enabled 
}

INITIAL {
	I_ou = 0
	active = -1
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
