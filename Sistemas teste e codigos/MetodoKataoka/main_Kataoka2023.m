clc; clear;

%% ---- Case Import ---- %%
%%%%Sistema de teste (CARTAO NO FORMATO ADOTADO PELO MATHEUS!)
%casobase = case14;
%casobase = case30christy';
%casobase = case118
%[NB,NA,~,pq,pv,ref,gen,Ybus,Vm,Va,Pg0,Qg,Pd0,Qd0,Qmax,Qmin,Vmax,Vmin,Sb,names,buses] = data(casobase);

%%%%Sistema de teste (CARTAO NO FORMATO ORIGINAL DO IEEE)
%system = 'ieee14.cdf';
%system = 'ieee30.cdf';
%system = 'ieee57.cdf';
system = 'ieee118.cdf';
% system = 'ieee300cdf.txt';
[NB,NA,~,pq,pv,ref,gen,Ybus,Vm,Va,Pg0,Qg,Pd0,Qd0,Qmax,Qmin,Vmax,Vmin,Sb,names,buses] = systemData(system);

npv = length(pv); npq = length(pq);         %Vector Length
nbus = npv+npq+1;

%% ---- Initial Values ---- %%
t = 0;                                      %Load Growth Parameter
wP = zeros(size(Vm)); wQ = zeros(size(Vm));
wP([pv; pq]) = 1; wQ(pq) = 1;                %Initial Eigenvector
it = 0;       %Iteration Counter
tol = 1e-4;   %Tolerance (Newton-Raphson)
s = 1;        %Initial Convergence Error
it_max = 20;

%% Load uncertainty initialization
zl(:,1) = Pd0/norm(Pd0);
sigma = 0.1*norm(Pd0);
r = sigma;

fpart = Pg0/sum(Pg0);
tan0 = Qd0/Pd0;

for caseWorst = 1 : 1 : 2
    if caseWorst == 1
        fprintf('Solution for Minimum- Critical Points \n')
        i = 1;
    elseif caseWorst == 2
        fprintf('\nSolution for Minimum-Total-Load Critical Point \n')
        i = 1;
    end
    while true
        i = i + 1;
        while true
            Pd=(1+t).*Pd0+r*zl(:,i-1); %Demanded Load, r = sigma * t^(1/2)
            Qd=Qd0.*Pd; %Qd(Pd) = Qd0*pd
            Pg=Pg0+fpart*(sum(t.*Pd0+r*zl(:,i-1)));
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

            g = df_dx.'*[wP([pv; pq]); wQ(pq)];     %Function g

            %% ---- Function h(u) Building ---- %%
            h = sum([wP([pv; pq]); wQ(pq)].*[wP([pv; pq]); wQ(pq)])-1;    %Function h

            %% ---- Function F(x,t,u) Building ---- %%
            F = [f; g; h];                                                %Function F
            %         dPd_dt = -fpg*sum(Pd0+sigma*zl(:,i-1)/(2*sqrt(t)))+Pd0+sigma*zl(:,i-1)/(2*sqrt(t)); dQd_dt = Qd0.*dPd_dt;
            dPd_dt = Pd0-fpart*sum(Pd0); dQd_dt = Qd0.*dPd_dt;
            df_dt = [dPd_dt([pv;pq]); dQd_dt(pq)];                        %Jacobian Matrix df_dt
            df_du = zeros(size(df_dx,1),size([wP([pv; pq]); wQ(pq)],1));  %Jacobian Matrix df_du

            %% ---- Function g(x,t,u) = Dxf.'u derivatives ---- %%
            [GAA1, GAV1, GVA1, GVV1] = d2Sbus_dV2(Ybus, V, wP);           %Complex Jacobian Matrix Dxxf.*uP
            [GAA2, GAV2, GVA2, GVV2] = d2Sbus_dV2(Ybus, V, wQ);           %Complex Jacobian Matrix Dxxf.*uQ

            M1 = [GAA1([pv; pq], [pv; pq]), GAV1([pv; pq], pq); GVA1(pq, [pv; pq]), GVV1(pq, pq)];
            M2 = [GAA2([pv; pq], [pv; pq]), GAV2([pv; pq], pq); GVA2(pq, [pv; pq]), GVV2(pq, pq)];

            dg_dx = real(M1)+imag(M2);
            dg_dt = zeros(size(df_dt));
            dg_du = df_dx.';

            %% ---- Function h(u) = u*u.'- 1 derivatives ---- %%
            dh_dx = zeros(1,size(df_dx,1));
            dh_dt = 0;
            dh_du = 2*[wP([pv; pq]); wQ(pq)].';

            %% ---- NR Loop ---- %%
            J = [df_dx, df_dt, df_du; dg_dx, dg_dt, dg_du; dh_dx, dh_dt, dh_du];    %Complete Jacobian Matrix
            x = [Va([pv; pq]); Vm(pq); t; wP([pv; pq]); wQ(pq)];                    %Previous State Vector

            s = J\F;      %Convergence Error
            x = x-s;      %Next State Vector

            %% ---- State Variable Update ---- %%
            Va([pv; pq]) = x(1:npv+npq);
            Vm(pq) = x(npv+npq+1:npv+2*npq);
            t = x(npv+2*npq+1);
            wP([pv; pq]) = x(npv+2*npq+2:2*npv+3*npq+1);
            wQ(pq) = x(2*npv+3*npq+2:end);
            it = it+1;
            if (norm(s)<=tol)||(it>it_max)
                break;
            end
        end
        w = [wP([pv; pq]);wQ(pq)];
        r = sigma*sqrt(t);
        aux1 = diag(ones(size([pv; pq])));
        aux2 = diag(Qd0(pq));
        df_dpl = blkdiag(aux1,aux2);
        aux3 = (w.')*df_dpl;
        if caseWorst == 1 %
            zl(2:end,i) = aux3([pv;pq]);
            %Acrescentado no código depois /norm(aux3([pv;pq]) para ficar como em Jabr_2009 (doi:10.1016/j.ijepes.2009.01.009)
            %Resultado dessa normalização adiiconal: o valor de t fica ainda menor
            zl(2:end,i) = zl(2:end,i)/norm(aux3([pv;pq]));
        elseif caseWorst == 2 %normal vector to the constant-total load surface
            nt = Pd/norm(Pd);
            nt(ref) = [];
            nc = (aux3([pv;pq])/norm(aux3([pv;pq])))'; %normal vector to limit surface
            F = nt'*nc; %nt'*nc
            E = nt'*Pd0(2:end); %nt'*dPls_dt
            D = nc'*Pd0(2:end); %nc'*dPls_dt
            C = (1/2*t^(-1/2))^2-D^2; %(dr_dt)^2-D^2;
            B = 2*(1/2*t^(-1/2))*(E-D*F); %2*dr_dt*(E-D*F)
            A = D^2 + E^2 - 2*D*E*F;
            syms a b zl2;
            b = min(double(solve(b^2*A+b*B+C,b)));
            a = double(solve((a*nc+b*nt)'*(Pd0(2:end))==(-1/2*t^(-1/2)),a));
            zl(2:end,i) = a*nc + b*nt;           
        end           
        t_values(i) = t; %t values
        %Total load power is expressed without adding the uncertainty
        %component (Pdo*(1+t))
        fprintf('Iteration #%d: zl norm = %f, t value = %.2f, and total load power (Pdo*(1+t)) = %.2f MW \n', i-1, norm(zl(:,i)), t,sum(Pd0*(1+t))*Sb)
        if norm(zl(:,i) - zl(:,i-1)) < 1e-8
            break;
        end
    end
end
