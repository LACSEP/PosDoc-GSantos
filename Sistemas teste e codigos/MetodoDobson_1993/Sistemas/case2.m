function mpc = sist3

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 1;

%% bus data
%	bus_i	type	Pd	Qd	 Gs	Bs	area   Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	 1	      3	     0	 0	 0	0	 1	    1	0	230	      1	     1.1	0.9;
	 2	      1  	0.5	0.3	 0	0	 1	    0.908	-7.9068	230	      1	     1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	 0	 0	 0	     0	    1	100	      1	     0	     0	0	0	0	0	0	0	0	0	0	0	0;	
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	  1	     2	  0.00	0.25	0.00	0	0	0	0	0	1	-360	360;	  
];
