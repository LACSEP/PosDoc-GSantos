clc; clear;

%% ---- Case Import ---- %%
%%%%Sistema de teste (CARTAO NO FORMATO ADOTADO PELO MATHEUS!)
tic;
addpath('.\Sistemas')

%%%%Sistema de teste (CARTAO NO FORMATO ORIGINAL DO IEEE)
%system = 'ieee14.cdf';
%system = 'ieee30.cdf';
%system = 'ieee57.cdf';
%system = '107barra.cdf';
%system = 'ieee118.cdf';
%system = 'NETS_NYPS_vb.cdf';
system = '107barras_ufes.cdf';
[NB,NA,~,pq,pv,ref,gen,Ybus,Vm,Va,Pg0,Qg,Pd0,Qd0,Qmax,Qmin,Vmax,Vmin,Sb,names,buses] = systemData(system);

%% ---- Ordering and size of data vectors ---- %%
npv = length(pv); npq = length(pq);
nbus = npv + npq + 1;
gen = sort(gen); nger = length(gen);

%% ---- Initial Values ---- %%
t = 0; %Load Growth Parameter
wP = ones(nbus,1)/2;
wP(ref) = 0;
wQ = ones(nbus,1);
wG = ones(nger,1);
Vr = Vm;
Qspc = Qg(gen);
if Qmax(ref) == 0
    Qmax(ref) = 999.99;
end
if Qmin(ref) == 0
    Qmin(ref) = -999.99;
end

%% ---- Control Constants ---- %%
sigk = 1e8; sigv = 1e-6; sigq = 1e-6;
tepa = 1e-8; tepr = 1e-8; qlst = 0.04;

%% ---- Matrices without changes ---- %%
dP_dQg = zeros(nbus-1,nger);
dQ_dQg = zeros(nbus,nger);
for i = 1 : 1 : nger
    j = gen(i);
    dQ_dQg(j,i) = -1;
end
dy_dVa = zeros(nger,nbus-1);
dA_dQg = zeros(nbus-1,nger);
dC_dVa = zeros(nger,nbus-1);

itFP = 0;
while true
    itFP = itFP + 1;
    %% ---- Smoothed reactive control ---- %%
    [diffy_Qg,diffy_v,diffy_Qg_Qg,diffy_Qg_v,diffy_v_v,diffy_v_Qg,Y] = CtrlQlims(sigk,sigv,sigq,Vm,Vr,Qspc,Qmax,Qmin,gen,nger);
    dy_dVm = zeros(nger,nbus);
    dy_dQg = zeros(nger, nger);
    for i = 1 : 1 : nger
        j = gen(i);
        dy_dVm(i,j) = diffy_v(i);
        % if Qg(j) > (Qmax(j) - sigq) || Qg(j) < (Qmin(j) + sigq)
        dy_dQg(i,i) = diffy_Qg(i);
        %  end
    end
    %% ---- Calculation of specified and calculated power ---- %%
    Qg(gen) = Qspc;
    Sesp = (Pg0-Pd0)+1i*(Qg-Qd0);          %Specified complex power
    V = Vm.*exp(1i*Va);                 %Voltage phasor
    I = Ybus*V;                         %Current phasor
    Scalc = diag(V)*conj(I);            %Calculated complex power
    %% ---- Jacobian function f(x,t,u) and its partial derives ---- %%
    [dS_dVa, dS_dVm] = dSbus_dV(Ybus, V);   %Complex derivatives
    A11 = real(dS_dVa); A11(ref,:) = []; A11(:,ref) = []; %dP_dVa
    A12 = real(dS_dVm); A12(ref,:) = [];    %dP_dVm
    A21 = imag(dS_dVa); A21(:,ref) = [];    %dQ_dVa
    A22 = imag(dS_dVm);                     %dQ_dVm
    df_dx = [A11 A12 dP_dQg; A21 A22 dQ_dQg; dy_dVa dy_dVm dy_dQg]; %df_dx matrix
    %% ---- Calculation of residues ---- %%
    eP = real(Scalc-Sesp);              %Active power mismatch
    eQ = imag(Scalc-Sesp);              %Reactive power mismatch
    ey = Y;                             %Control mismatch
    if norm(eP) < tepa || norm(eQ) < tepr
        break;
    end
    f = [eP(1:end ~= ref);eQ;ey];  %Function f
    s = mldivide(df_dx,f);      %Convergence Error
    x = [Va(1:end ~= ref); Vm; Qspc];
    x = x-s;      %Next State Vector
    Va(1:end ~= ref) = x(1:nbus-1);
    Vm = x(nbus:2*nbus-1);
    Qspc = x(2*nbus:2*nbus+npv);
end
Pg0 = real(Scalc)+Pd0;

tol = 1e-7;   %Tolerance (Newton-Raphson)
itMax = 200;
it = 0;       %Iteration Counter
kpq = zeros(nbus,1); kpq(pq) = 1;

while true
    %% ---- Load increment ---- %%
    Pd = (1+t).*Pd0; %Demanded active load
    Qd = (1+t).*Qd0; %Demanded reactive load
    Pg = Pg0; %ATTENTION! The participation factor is not included!
    Qg(gen) = Qspc;
    %% ---- Smoothed reactive control ---- %%
    [diffy_Qg,diffy_v,diffy_Qg_Qg,diffy_Qg_v,diffy_v_v,diffy_v_Qg,Y] = CtrlQlims(sigk,sigv,sigq,Vm,Vr,Qspc,Qmax,Qmin,gen,nger);
    dy_dVm = zeros(nger,nbus);
    dy_dQg = zeros(nger, nger);
    for i = 1 : 1 : nger
        j = gen(i);
        dy_dVm(i,j) = diffy_v(i);
        % if Qg(j) > (Qmax(j) - sigq) || Qg(j) < (Qmin(j) + sigq)
        dy_dQg(i,i) = diffy_Qg(i);
        %  end
    end
    dy_dvdv = zeros(nbus,nbus);
    dB_dQg = zeros(nbus,nger);
    for i = 1 : 1 : nger
        j = gen(i);
        dy_dvdv(j,j) = -diffy_v_v(i)*wG(i);
        dB_dQg(j,i) = -diffy_v_Qg(i)*wG(i);
    end
    dC_dVm = zeros(nger,nbus);
    dC_dQg = zeros(nger,nger);
    for i = 1 : 1 : nger
        j = gen(i);
        %    if Qg(j) > (Qmax(j) - sigq) || Qg(j) < (Qmin(j) + sigq)
        dC_dVm(i,j) = -diffy_Qg_v(i)*wG(i);
        dC_dQg(i,i) = -diffy_Qg_Qg(i)*wG(i);
        %    end
    end

    %% ---- Calculation of residues ---- %%
    Sesp = (Pg-Pd)+1i*(Qg-Qd);          %Specified complex power
    V = Vm.*exp(1i*Va);                 %Voltage phasor
    I = Ybus*V;                         %Current phasor
    Scalc = diag(V)*conj(I);            %Calculated complex power
    eP = real(Scalc-Sesp);              %Active power mismatch
    eQ = imag(Scalc-Sesp);              %Reactive power mismatch
    ey = Y;                             %Control mismatch

    %% ---- Jacobian function f(x,t,u) and its partial derives ---- %%
    [dS_dVa, dS_dVm] = dSbus_dV(Ybus, V);   %Complex derivatives
    A11 = real(dS_dVa); A11(ref,:) = []; A11(:,ref) = []; %dP_dVa
    A12 = real(dS_dVm); A12(ref,:) = [];    %dP_dVm
    A21 = imag(dS_dVa); A21(:,ref) = [];    %dQ_dVa
    A22 = imag(dS_dVm);                     %dQ_dVm
    f = [eP(1:end ~= ref);eQ;ey];  %Function f
    df_dx = [A11 A12 dP_dQg; A21 A22 dQ_dQg; dy_dVa dy_dVm dy_dQg]; %df_dx matrix
    dPd_dt = Pd0; dQd_dt = Qd0;
    df_dt = [dPd_dt(1:end ~= ref); dQd_dt; zeros(nger,1)];  %df_dt matrix
    df_du = zeros(nbus*2+nger-1,nbus*2+nger-1); %df_du matrix

    %% ---- Hessian function g(x,t,u) and its partial derives ---- %%
    g = df_dx.' * [wP(1:end ~= ref); wQ; wG];    %Function g
    [GAA1, GAV1, GVA1, GVV1] = d2Sbus_dV2(Ybus, V, wP);           %Complex jacobian matrix Dxxf.*uP
    [GAA2, GAV2, GVA2, GVV2] = d2Sbus_dV2(Ybus, V, wQ);           %Complex jacobian matrix Dxxf.*uQ
    dA_dVa = real(GAA1) + imag(GAA2); dA_dVa(ref,:) = []; dA_dVa(:,ref) = [];
    dA_dVm = real(GAV1) + imag(GAV2); dA_dVm(ref,:) = [];
    dB_dVa = real(GVA1) + imag(GVA2); dB_dVa(:,ref) = [];
    dB_dVm = real(GVV1) + imag(GVV2) + dy_dvdv;
    dg_dx = [dA_dVa dA_dVm dA_dQg; dB_dVa dB_dVm dB_dQg; dC_dVa, dC_dVm, dC_dQg];
    dg_dt = zeros(size(df_dt));
    dg_du = df_dx.';

    %% ---- Function h(u) and its partial derives ---- %%
    h = sum([wP(1:end ~= ref); wQ; wG] .* [wP(1:end ~= ref); wQ; wG])-1;    %Function h
    dh_dx = zeros(1,size(df_dx,1));
    dh_dt = 0;
    dh_du = 2*[wP(1:end ~= ref); wQ; wG].';


    %% ---- Updating state variables ---- %%
    F = [f; g; h];  %Function F
    J = [df_dx, df_dt, df_du; dg_dx, dg_dt, dg_du; dh_dx, dh_dt, dh_du];    %Complete Jacobian Matrix
    JF = full(J);
    x = [Va(1:end ~= ref); Vm; Qspc; t; wP(1:end ~= ref); wQ; wG];                  %Previous State Vect
    s = mldivide(JF,F);      %Convergence Error
    x = x-s;      %Next State Vector
    Va(1:end ~= ref) = x(1:nbus-1);
    Vm = x(nbus:2*nbus-1);
    Qspc = x(2*nbus:2*nbus+npv);
    t = x(2*nbus+nger);
    wP(1:end ~= ref) = x(2*nbus+nger+1:3*nbus+npv);
    wQ = x(3*nbus+npv+1:4*nbus+npv);
    wG = x(4*nbus+npv+1:end);

    it = it+1;

    if (norm(s)<=tol)||(it>itMax)
        break;
    end
end
toc;