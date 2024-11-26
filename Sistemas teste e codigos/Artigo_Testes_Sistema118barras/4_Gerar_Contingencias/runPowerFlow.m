function [Vm, Va, Pg0, Qg, itPF] = runPowerFlow(Vm,Va,Vr,Pd0,Qd0,Pg0,Qg,nbus,npv,ref,gen,nger,Qspc,Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,Ybus,tepa,tepr,itPFmax)
itPF = 0;  
while true
    itPF = itPF + 1;
    %% ---- Smoothed reactive control ---- %%
    [diffy_Qg,diffy_v,Y] = CtrlQlimsPF(1e8,1e-6,1e-6,Vm,Vr,Qspc,Qmax,Qmin,gen,nger);
    dy_dVm = zeros(nger,nbus);
    dy_dQg = zeros(nger, nger);
    for i = 1 : 1 : nger
        j = gen(i);
        dy_dVm(i,j) = diffy_v(i);
        dy_dQg(i,i) = diffy_Qg(i);
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
    if itPF > itPFmax
      %  fprintf('Fluxo de potência não convergiu com %d iteracoes', itPFmax);
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
Pg0(Pg0<10e-4) = 0;
Qg(Qg<10e-4) = 0;