clc; clear;

warning off;

addpath('.\Sistemas');

%% ---- Case Import ---- %%
%%%%Sistema de teste (CARTAO NO FORMATO ADOTADO PELO MATHEUS!)
% casobase = case2;
 casobase = case5;
% casobase = case14;
% casobase = case30christy';
[NB,NA,~,pq,pv,ref,gen,Ybus,Vm,Va,Pg0,Qg,Pd0,Qd0,Qmax,Qmin,Vmax,Vmin,Sb,names,buses] = data(casobase);

npv = length(pv); npq = length(pq);         %Vector Length
nbus = npv+npq+1;

%% ---- Inicialization ---- %%

it = 0;       %Iteration Counter
tol = 1e-4;   %Tolerance (Newton-Raphson)
s = 1;        %Initial Convergence Error
it_max = 30;

wP = zeros(size(Vm)); wQ = zeros(size(Vm));
wP([pv; pq]) = 1; wQ(pq) = 1;                %Initial Eigenvector
t = 1;

idx = 1;
while true
    if idx == 1               
       n0 = [Pd0/norm(Pd0); Qd0/norm(Qd0)];
    end
    while true
        it = it + 1;
        %% ---- Function f(x,t) Building ---- %%
        Pd = Pd0+t*n0(1:nbus,1); Qd = Qd0+t*n0(nbus+1:end,1);        
        Pg = Pg0;
        Sesp = (Pg-Pd)+1i*(Qg-Qd);          %Specified Complex Power
        V = Vm.*exp(1i*Va);                 %Voltage Phasor
        I = Ybus*V;                         %Current Phasor
        Scalc = diag(V)*conj(I);            %Calculated Complex Power
        eP = real(Scalc-Sesp);              %Active Power Mismatch
        eQ = imag(Scalc-Sesp);              %Reactive Power Mismatch
        f = [eP([pv; pq]);eQ(pq)];          %Function f

        %% ---- Function g(x,t,u) Building ---- %%
        [dS_dVa, dS_dVm] = dSbus_dV(Ybus, V);   %Complex Derivatives
        A11 = real(dS_dVa([pv; pq], [pv; pq])); %dP_dVa
        A12 = real(dS_dVm([pv; pq], pq));       %dP_dVm
        A21 = imag(dS_dVa(pq, [pv; pq]));       %dQ_dVa
        A22 = imag(dS_dVm(pq, pq));             %dQ_dVm
        df_dx = [ A11 A12 ; A21 A22 ];          %Jacobian Matrix df_dx
        w = [wP([pv; pq]);wQ(pq)];
        g = df_dx.'*w;     %Function g

        %% ---- Function h(w) Building ---- %%
        h = sum([w].*[w])-1;    %Function h

        %% ---- Function F(x,t,w) Building ---- %%
        F = [f; g; h];                                                %Function F
        %         dPd_dt = -fpg*sum(Pd0+sigma*zl(:,i-1)/(2*sqrt(t)))+Pd0+sigma*zl(:,i-1)/(2*sqrt(t)); dQd_dt = Qd0.*dPd_dt;
        dPd_dt = n0(1:nbus,1); dQd_dt = n0(nbus+1:end,1);
        df_dt = [dPd_dt([pv;pq]); dQd_dt(pq)];                        %Jacobian Matrix df_dt
        df_du = zeros(size(df_dx,1),size(w,1));  %Jacobian Matrix df_du

        %% ---- Function g(x,t,w) = Dxf.'w derivatives ---- %%
        [GAA1, GAV1, GVA1, GVV1] = d2Sbus_dV2(Ybus, V, wP);           %Complex Jacobian Matrix Dxxf.*wP
        [GAA2, GAV2, GVA2, GVV2] = d2Sbus_dV2(Ybus, V, wQ);           %Complex Jacobian Matrix Dxxf.*wQ

        M1 = [GAA1([pv; pq], [pv; pq]), GAV1([pv; pq], pq); GVA1(pq, [pv; pq]), GVV1(pq, pq)];
        M2 = [GAA2([pv; pq], [pv; pq]), GAV2([pv; pq], pq); GVA2(pq, [pv; pq]), GVV2(pq, pq)];

        dg_dx = real(M1)+imag(M2);
        dg_dt = zeros(size(df_dt));
        dg_dw = df_dx.';

        %% ---- Function h(u) = u*u.'- 1 derivatives ---- %%
        dh_dx = zeros(1,size(df_dx,1));
        dh_dt = 0;
        dh_dw = 2*w.';

        %% ---- NR Loop ---- %%
        J = [df_dx, df_dt, df_du; dg_dx, dg_dt, dg_dw; dh_dx, dh_dt, dh_dw];    %Complete Jacobian Matrix
        x = [Va([pv; pq]); Vm(pq); t; wP([pv; pq]); wQ(pq)];                    %Previous State Vector

        s = J\F;      %Convergence Error
        x = x-s;      %Next State Vector

        %% ---- State Variable Update ---- %%
        Va([pv; pq]) = x(1:npv+npq);
        Vm(pq) = x(npv+npq+1:npv+2*npq);
        t = x(npv+2*npq+1);
        wP([pv; pq]) = x(npv+2*npq+2:2*npv+3*npq+1);
        wQ(pq) = x(2*npv+3*npq+2:end);
        if (norm(s)<=tol)||(it>it_max)
            break;
        end
    end    
    n0v(:,idx) = n0;
    idx = idx + 1;
    n0([pv;pq]) = wP([pv;pq]);
    n0(nbus+pq) = wQ(pq);
    n0v(:,idx) = n0;
    fprintf('Iteracao %d \n', idx-1);    
    lambP(:,idx) = round(t*n0(1:nbus,1),3);
    idx_minP = 0;
    value_minP = 0;
    for i = 1 : 1 : npv
        if idx_minP == 0
            idx_minP = pv(i);
            value_minP = lambP(pv(i),idx);
        else
            if value_minP > lambP(pv(i),idx)
                idx_minP = pv(i);
                value_minP = lambP(pq(i),idx);
            end
        end
    end
    for i = 1 : 1 : npq
        if idx_minP == 0
            idx_minP = pq(i);
            value_minP = lambP(pq(i),idx);
        else
            if value_minP > lambP(pq(i),idx)
                idx_minP = pq(i);
                value_minP = lambP(pq(i),idx);
            end
        end
    end
    fprintf('Mínimo P-P0 na barra %d com valor igual a %f MW \n', idx_minP, value_minP*Sb);
  
    lambQ(:,idx) = round(t*n0(nbus+1:end,1),3);
    idx_minQ = 0;
    value_minQ = 0;    
    for i = 1 : 1 : npq
        if idx_minQ == 0
            idx_minQ = pq(i);
            value_minQ = lambQ(pq(i),idx);
        else
            if value_minQ > lambQ(pq(i),idx)
                idx_minQ = pq(i);
                value_minQ = lambQ(pq(i),idx);
            end
        end
    end    
    fprintf('Mínimo Q-Q0 na barra %d com valor igual a %f MVAr \n', idx_minQ, value_minQ*Sb);
    if norm(n0v(:,idx)-n0v(:,idx-1))<1e-4            
        break;
    end    
    
end



