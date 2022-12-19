NEURON {
	POINT_PROCESS RxDSyn
	RANGE tau
}

PARAMETER {
	tau = 0.1 (ms) <1e-9,1e9>
}

STATE { g }

INITIAL { g=0 }

BREAKPOINT {SOLVE state METHOD cnexp}

DERIVATIVE state { g' = -g/tau }
NET_RECEIVE(weight) {
	g = g + weight
}

