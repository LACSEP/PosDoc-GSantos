function [t, dPtotal, it] = runDirectMethod(Vm,Va,Vr,Pd0,Pg0,Qg,nbus,npv,ref,gen,nger,Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,dA_dQg,dC_dVa,Ybus,bpl,tan0,fpart,tol,itPoCmax)                            
%% ---- Initial Values ---- %%
t = 0;
wP = ones(nbus,1);
wP(ref) = 0;
wQ = ones(nbus,1);
wG = ones(nger,1);

bpg = fpart*sum(bpl);

it = 0;
while true
    %% ---- Load increment ---- %%
    Pd = Pd0 + t*bpl; %Demanded active load
    Qd = Pd.*tan0;
    Pg = Pg0 + t*bpg; 
    %% ---- Smoothed reactive control ---- %%
    [diffy_Qg,diffy_v,diffy_Qg_Qg,diffy_Qg_v,diffy_v_v,diffy_v_Qg,Y] = CtrlQlims(1e8,1e-6,1e-6,Vm,Vr,Qg(gen),Qmax,Qmin,gen,nger);
    dy_dVm = zeros(nger,nbus);
    dy_dQg = zeros(nger, nger);
    for i = 1 : 1 : nger
        j = gen(i);
        dy_dVm(i,j) = diffy_v(i);
        dy_dQg(i,i) = diffy_Qg(i);
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
        dC_dVm(i,j) = -diffy_Qg_v(i)*wG(i);
        dC_dQg(i,i) = -diffy_Qg_Qg(i)*wG(i);
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
    dPd_dt = bpl-bpg; dQd_dt = tan0.*bpl;
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
    x = [Va(1:end ~= ref); Vm; Qg(gen); t; wP(1:end ~= ref); wQ; wG];                  %Previous State Vect
    s = mldivide(JF,F);      %Convergence Error
    
    x = x-s;      %Next State Vector
    Va(1:end ~= ref) = x(1:nbus-1);
    Vm = x(nbus:2*nbus-1);
    Qg(gen) = x(2*nbus:2*nbus+npv);
    t = x(2*nbus+nger);
    wP(1:end ~= ref) = x(2*nbus+nger+1:3*nbus+npv);
    wQ = x(3*nbus+npv+1:4*nbus+npv);
    wG = x(4*nbus+npv+1:end);

    it = it+1;  

    if (norm(s)<=tol)||(it>itPoCmax)
        if it > itPoCmax
       %     fprintf('PoC não convergiu \n');
        end
        break;
    end
end
dPtotal = sum(bpl*t);
%toc;