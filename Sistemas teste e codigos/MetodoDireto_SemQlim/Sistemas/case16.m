function mpc = case16
%CASE16 Power sistema-teste brasileiro reduzido com 16 barras

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	0	0	0	0	1	0.97	0.52	0	1	1.06	0.94;
	2	1	0	0	0	0	1	0.99	-1.30	0	1	1.06	0.94;
	3	1	26	16	0	-90	1	1.00	-4.20	0	1	1.06	0.94;
	4	1	58	39	0	0	2	1.00	-5.60	0	1	1.06	0.94;
	5	1	16	11	0	0	2	1.01	-4.80	0	1	1.06	0.94;
	6	1	10	7	0	0	2	1.00	-3.20	0	1	1.06	0.94;
	7	1	9	8	0	0	2	1.01	-4.90	0	1	1.06	0.94;
	8	1	21	13	0	0	2	1.02	-6.00	0	1	1.06	0.94;
	9	2	0	0	0	0	2	0.97	-6.00	0	1	1.06	0.94;
	10	1	23	16	0	20 	1	0.97	-6.80	0	1	1.06	0.94;
	11	1	22	15	0	-30	1	0.99	-6.70	0	1	1.06	0.94;
	12	1	56	31	0	0	1	1.01	-6.00	0	1	1.06	0.94;
	13	1	33	19	0	-30	1	1.01	-6.30	0	1	1.06	0.94;
	14	1	32	21	0	0	1	1.01	-5.70	0	1	1.06	0.94;
    15	1	0	0	0	0	1	1.01	-4.80	0	1	1.06	0.94;
    16	3	0	0	0	0	2	0.97	0.00	0	1	1.06	0.94;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	200 -130  180 -180	  1.06	100	1	332.4	0	0	0	0	0	0	0	0	0	0	0	0;
	9	0	-39	  70  -50	  1.045	 100	1	140	0	0	0	0	0	0	0	0	0	0	0	0;
	16	108	-53	  180 -180	  1.01	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0.0000	0.0150	0.0000	0	0	0	0	0	1	-360	360;
	2	3	0.0023	0.0248	1.6960	0	0	0	0	0	1	-360	360;
	3	12	0.0024	0.0264	0.4505	0	0	0	0	0	1	-360	360;
	3	14	0.0000	0.0833	0.0000	0	0	0	0.9742	0	1	-360	360;
	3	15	0.0040	0.0419	0.6126	0	0	0	0	0	1	-360	360;
	4	5	0.0034	0.0368	0.5390	0	0	0	0	0	1	-360	360;
	4 	15	0.0000	0.0667	0.0000	0	0	0	0.9991	0	1	-360	360;
	5	6	0.0050	0.0536	0.1960	0	0	0	0	0	1	-360	360;
	6	7	0.0060	0.0637	0.2328	0	0	0	0	0	1	-360	360;
	6	16	0	0.5000	0	0	0	0	0.932	0	1	-360	360;
	7	8	0.0047	0.0503	0.1838	0	0	0	0	0	1	-360	360;
	8	9	0.0000	0.1250	0.0000	0	0	0	0	0	1	-360	360;
	10	8	0.0000	0.1000	0.0000	0	0	0	0.9030	0	1	-360	360;
	10	11	0.0034	0.0372	0.6360	0	0	0	0	0	1	-360	360;
	11	12	0.0039	0.0434	0.7420	0	0	0	0	0	1	-360	360;
	12	13	0.0011	0.0124	0.8480	0	0	0	0	0	1	-360	360;
];

%% bus names
mpc.bus_name = {
	'GERADOR-1';
	'BARRA-2';
	'BARRA-3';
	'BARRA-4';
	'BARRA-5';
	'BARRA-6';
	'BARRA-7';
	'BARRA-8';
	'C.SINCRONO';
	'BARRA-10';
	'BARRA-11';
	'BARRA-12';
	'BARRA-13';
	'BARRA-14';
    'BARRA-15';
    'GERADOR-2';
};