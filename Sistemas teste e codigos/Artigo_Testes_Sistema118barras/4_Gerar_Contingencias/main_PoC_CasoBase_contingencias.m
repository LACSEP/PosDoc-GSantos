clc; clear;

warning off;

%% Dados do sistema teste
contigencias = [1	2	3	4	5	6	7	8	10	11	12	13	14	15	16	17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32	33	34	35	36	37	38	39	40	41	42	43	44	45	46	47	48	49	50	51	52	53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68	69	70	71	72	73	74	75	76	77	78	79	80	81	82	83	84	85	86	87	88	89	90	91	92	93	94	95	96	97	98	99	100	101	102	103	104	105	106	107	108	109	110	111	112	114	115	116	117	118	119	120	121	122	123	124	125	126	127	128	129	130	131	132	133	135	136	137	138	139	140	141	142	143	144	145	146	147	148	149	150	151	152	153	154	155	156	157	158	159	160	161	162	163	164	165	166	167	168	169	170	171	172	173	174	175	178	179	180	181	182	184	185	186];
nContigencias = length(contigencias);

addpath('.\Sistemas')

system = 'case118bus_modificado.cdf';
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,dbranch] = systemData(system);

system = [];

%Ramos do sistema
nBranch = length(dbranch);

tic;
for caso = 1 : 1 : nContigencias
    %% Sistema teste com contigências
    system = strcat('Linha',num2str(contigencias(caso)),'_',num2str(dbranch(contigencias(caso)).iBus),'_',num2str(dbranch(contigencias(caso)).jBus),'.cdf');
    [NB,NA,~,pq,pv,ref,gen,Ybus,Vm,Va,Pg0,Qg,Pd0,Qd0,Qmax,Qmin,Vmax,Vmin,Sb,names,buses,~] = systemData(system);

    %Tensão de referência
    Vr = Vm;

    %Potência gerada de referência
    Pg0_old = Pg0;
    Pgt_old = sum(Pg0)-Pg0(ref);

    %Barras pv, pq e geradores (inclui slack)
    npv = length(pv); npq = length(pq);
    gen = sort(gen); nger = length(gen);
    nbus = npv + npq + 1;

    %Limite de reativo da barra de referência
    if Qmax(ref) == 0
        Qmax(ref) = 999.99;
    end
    if Qmin(ref) == 0
        Qmin(ref) = -999.99;
    end

    %Eolicas
    wind = find(Pd0<0);  nWind = length(wind); %Buses with wind farm

    %Barras de carga
    loadBuses = []; j = 0;
    for i = 1 : 1 : nbus
        if Pd0(i) > 0
            j = j + 1;
            loadBuses = [loadBuses; buses(i)];
        end
    end
    nLoad = length(loadBuses); %Buses with load
    tan0 = zeros(nbus,1); tan0(loadBuses) = Qd0(loadBuses)./Pd0(loadBuses);
    Pd0_total = sum(Pd0(loadBuses));

    %Fator de participação dos geradores
    fpart = zeros(nbus,1);
    fpart(Pg0>0) = 0.047874377633091;
    fpart(ref) = 0.138261202604366;

    %% Control matrices without changes

    dP_dQg = zeros(nbus-1,nger);
    dQ_dQg = zeros(nbus,nger);
    for i = 1 : 1 : nger
        j = gen(i);
        dQ_dQg(j,i) = -1;
    end
    dy_dVa = zeros(nger,nbus-1);
    dA_dQg = zeros(nbus-1,nger);
    dC_dVa = zeros(nger,nbus-1);

    %% Caso base

    % Power flow
    itPFmax = 100; tepa = 1E-6; tepr = 1E-6;
    [Vm, Va, Pg0, Qg, itPF] = runPowerFlow(Vm,Va,Vr,Pd0,Qd0,Pg0,Qg,nbus,npv,ref,gen,nger,Qg(gen),Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,Ybus,tepa,tepr,itPFmax);
    Pg0(Pg0<10E-4) = 0;
    Qg(Qg<10E-4) = 0;
    % Direct Method
    itPoCmax = 200; tol = 1E-6;
    bpl = Pd0; bpl(wind) = 0;
    Pd0total = sum(bpl);
    [tload, dPtotal, itPoC] = runDirectMethod(Vm,Va,Vr,Pd0,Pg0,Qg,nbus,npv,ref,gen,nger,Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,dA_dQg,dC_dVa,Ybus,bpl,tan0,fpart,tol,itPoCmax);

    fprintf('----------------------CONTINGENCIA NA LINHA %d----------------------\n',contigencias(caso))
    fprintf('Resultado do método direto para o caso base (sem incertezas): \n');
    fprintf('Convergência em %d iteracoes com acrescimo de carga igual a %.2f MW (%.2f%%), lambda = %.4f \n', itPoC, dPtotal*Sb, dPtotal/Pd0total*100, tload);

    resultados(caso).Linha = contigencias(caso);
    resultados(caso).BarraDe = dbranch(contigencias(caso)).iBus;
    resultados(caso).BarraPara = dbranch(contigencias(caso)).jBus;
    resultados(caso).IncrementoMW = dPtotal*Sb;
    resultados(caso).Incrementoper = dPtotal/Pd0total*100;
    resultados(caso).lambda = tload;
    resultados(caso).itPoC = itPoC;

end
