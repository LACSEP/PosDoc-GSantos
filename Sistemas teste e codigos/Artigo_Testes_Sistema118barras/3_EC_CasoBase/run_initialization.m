%% Dados do sistema teste

[NB,NA,area,pq,pv,ref,gen,Ybus,Vm,Va,Pg0,Qg0,Pd0,Qd0,Qmax,Qmin,Vmax,Vmin,Sb,names,buses,~] = systemData(system);

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
[Vm, Va, Pg, Qg, itPF] = runPowerFlow(Vm,Va,Vr,Pd0,Qd0,Pg0,Qg0,nbus,npv,ref,gen,nger,Qg0(gen),Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,Ybus,tepa,tepr,itPFmax);
Pg(Pg0<10E-4) = 0;
Qg(Qg<10E-4) = 0;
% Direct Method
itPoCmax = 200; tol = 1E-7;
bpl = Pd0; bpl(wind) = 0;
Pd0total = sum(bpl);
[tloadCB, dPtotal, itPoC] = runDirectMethod(Vm,Va,Vr,Pd0,Pg,Qg,nbus,npv,ref,gen,nger,Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,dA_dQg,dC_dVa,Ybus,bpl,tan0,fpart,tol,itPoCmax);

fprintf('----------------------CASO BASE----------------------\n')
fprintf('Resultado do método direto para o caso base (sem incertezas): \n');
fprintf('Convergência em %d iteracoes com acrescimo de carga igual a %.2f MW (%.2f%%), lambda = %.4f \n', itPoC, dPtotal*Sb, dPtotal/Pd0total*100, tloadCB);



%% Matriz de correlacao das cargas e eólica

%Cargas
muLoad = bpl(loadBuses)';
for i = 1 : 1 : nLoad
    for j = 1 : 1 : nLoad
        if i == j
            SigmaLoad(i,j) = 1;
        elseif area(loadBuses(i)) == area(loadBuses(j))
            SigmaLoad(i,j) = 0.8;
        elseif (area(loadBuses(i)) == 1 && area(loadBuses(j)) == 2) || (area(loadBuses(i)) == 2 && area(loadBuses(j)) == 1)
            SigmaLoad(i,j) = 0.6;
        elseif (area(loadBuses(i)) == 2 && area(loadBuses(j)) == 3) || (area(loadBuses(i)) == 3 && area(loadBuses(j)) == 2)
            SigmaLoad(i,j) = 0.6;
        elseif (area(loadBuses(i)) == 3 && area(loadBuses(j)) == 4) || (area(loadBuses(i)) == 4 && area(loadBuses(j)) == 3)
            SigmaLoad(i,j) = 0.6;
        elseif (area(loadBuses(i)) == 4 && area(loadBuses(j)) == 5) || (area(loadBuses(i)) == 5 && area(loadBuses(j)) == 4)
            SigmaLoad(i,j) = 0.6;
        elseif (area(loadBuses(i)) == 4 && area(loadBuses(j)) == 6) || (area(loadBuses(i)) == 6 && area(loadBuses(j)) == 4)
            SigmaLoad(i,j) = 0.6;
        elseif (area(loadBuses(i)) == 5 && area(loadBuses(j)) == 6) || (area(loadBuses(i)) == 6 && area(loadBuses(j)) == 5)
            SigmaLoad(i,j) = 0.6;
        else
            SigmaLoad(i,j) = 0.4;
        end
    end
end
SigmaLoad_pot = zeros(nLoad,nLoad);
sigLoad = 0.1*muLoad;
for i = 1 : 1 : nLoad
    for j = 1 : 1 : nLoad
        SigmaLoad_pot(i,j) =  SigmaLoad(i,j) * (0.1*Pd0(loadBuses(i))) * (0.1*Pd0(loadBuses(j)));
    end
end
%windBus = [8 32 42 55 76 92 105];
%Eolicas
SigmaWind = eye(nWind,nWind);
SigmaWind(1,2) = 0.3; SigmaWind(2,1) = 0.3; %Correlação entre as eólicas das barras 8 e 32
SigmaWind(3,4) = 0.3; SigmaWind(4,3) = 0.3; %Correlação entre as eólicas das barras 42 e 55
SigmaWind(3,5) = 0.3; SigmaWind(5,3) = 0.3; %Correlação entre as eólicas das barras 42 e 76
SigmaWind(4,5) = 0.3; SigmaWind(5,4) = 0.3; %Correlação entre as eólicas das barras 55 e 76
SigmaWind(6,7) = 0.3; SigmaWind(7,6) = 0.3; %Correlação entre as eólicas das barras 92 e 105
SigmaWind_pot = zeros(nWind,nWind);
muWind = 2*ones(1,nWind);
dp = 0.25; % 0.25*2= 0.5 pu (50 MW)
sigWind = 0.25*muWind;
for i = 1 : 1 : nWind
    for j = 1 : 1 : nWind
        SigmaWind_pot(i,j) = SigmaWind(i,j)*(dp*muWind(i))*(dp*muWind(j)); % 25% de desvio padrão
    end
end

mu = [muLoad muWind]; %Média das potências consumidas pelas cargas e geradas pelas eólicas 
sig = [sigLoad sigWind]; %Desvio padrão das potências consumidas pelas cargas e geradas pelas eólicas 
Sigma = [SigmaLoad_pot zeros(nLoad, nWind); zeros(nWind,nLoad) SigmaWind_pot]; %Matriz Sigma com os valores absolutos
corrSigma = [SigmaLoad zeros(nLoad, nWind); zeros(nWind,nLoad) SigmaWind]; %Matriz de correlação

nLW = length(mu);
nW1 = nLoad + 1;

