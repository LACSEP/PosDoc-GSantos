clc; clear;

warning off;

%% ---- Power system data and parameter initialization ---- %%
% casebase -> CARTAO NO FORMATO ADOTADO PELO MATHEUS!
casobase = case5;
[nbus,narea,area,pq,pv,ref,gen,Ybus,Vm0,Va0,Pg0,Qg0,Pd0,Qd0,Qmax2,Qmin2,Vmax,Vmin,Sb,names,buses] = data(casobase);

%Ordering and size of data vectors
npv = length(pv); npq = length(pq);
gen = sort(gen); nger = length(gen);

%Reactive power limit of SGs
Qmax = zeros(nbus,1); Qmax(gen) = Qmax2;
Qmin = zeros(nbus,1); Qmin(gen) = Qmin2;
if Qmax(ref) == 0   Qmax(ref) = 999.99;     end
if Qmin(ref) == 0   Qmin(ref) = -999.99;    end

%Initial Values
Vr = Vm0;

%% Control Constants
sigk = 1e8; sigv = 1e-6; sigq = 1e-6;
tepa = 1e-8; tepr = 1e-8; qlst = 0.04;

%% Matrices without changes
dP_dQg = zeros(nbus-1,nger);
dQ_dQg = zeros(nbus,nger);
for i = 1 : 1 : nger
    j = gen(i);
    dQ_dQg(j,i) = -1;
end
dy_dVa = zeros(nger,nbus-1);
dA_dQg = zeros(nbus-1,nger);
dC_dVa = zeros(nger,nbus-1);

%Initial Values
t = 0; %Load Growth Parameter
wP = ones(nbus,1); wP(ref) = 0;
wQ = ones(nbus,1);
wG = ones(nger,1);
tol = 1e-7;  itMax = 200; %Tolerance (tolerance and maximum number of iterations)

%% Caso base

% Power flow
itPFmax = 100; tepa = 1E-6; tepr = 1E-6;
[Vm, Va, Pg0, Qg, itPF] = runPowerFlow(Vm0,Va0,Vr,Pd0,Qd0,Pg0,Qg0,nbus,npv,ref,gen,nger,Qg0(gen),Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,Ybus,tepa,tepr,itPFmax);
Pg0(Pg0<10E-4) = 0;
Qg(Qg<10E-4) = 0;

%Dispatch of SGs
fpart = zeros(nbus,1);
fpart(gen) = Pg0(gen)/sum(Pg0);

%Load growth factor
bpl = Pd0;
loadBuses = find(Pd0>0); nLoad = length(loadBuses); %Loaded buses
tan0 = zeros(nbus,1); tan0(loadBuses) = Qd0(loadBuses)./Pd0(loadBuses);

%Initial Values
t = 0; %Load Growth Parameter
wP = ones(nbus,1); wP(ref) = 0;
wQ = ones(nbus,1);
wG = ones(nger,1);
tol = 1e-7;  itMax = 200; itPoCmax = 100; %Tolerance (tolerance and maximum number of iterations)

Pd0total = sum(bpl);
[tloadCB, dPtotal, itPoC] = runDirectMethod(Vm,Va,Vr,Pd0,Pg0,Qg,nbus,npv,ref,gen,nger,Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,dA_dQg,dC_dVa,Ybus,bpl,tan0,fpart,tol,itPoCmax);

fprintf('----------------------CASO BASE----------------------\n')
fprintf('Resultado do método direto para o caso base (sem incertezas): \n');
fprintf('Convergência em %d iteracoes com acrescimo de carga igual a %.2f MW (%.2f%%), lambda = %.4f \n', itPoC, dPtotal*Sb, dPtotal/Pd0total*100, tloadCB);

%% ---- Test case determination and MCS/CE method ---- %%
percViol = 0.07; iteste = 0;
for cteste = 1 : 1 : 3
    fprintf('---- TESTE %d ---- \n',cteste)
    if cteste == 1
        dp = 0.3; correlacao = 0;
    elseif cteste == 2
        dp = 0.3; correlacao = 1; rho = 0.4;
    elseif cteste == 3
        dp = 0.3; correlacao = 1; rho = 0.8;
    end

    %Load power average vector
    mu = Pd0(loadBuses)';
    sig = dp*Pd0(loadBuses);

    %Correlation matrix
    corrSigma = zeros(nLoad,nLoad); %Normalized values
    Sigma = zeros(nLoad,nLoad); %Absolute values
    for i = 1 : 1 : nLoad
        for j = 1 : 1 : nLoad
            if i == j
                corrSigma(i,j) = 1;
                Sigma(i,j) = dp*Pd0(loadBuses(i)) * dp*Pd0(loadBuses(j));
            else
                if correlacao ~= 0
                    corrSigma(i,j) = rho;
                    Sigma(i,j) = rho*dp*Pd0(loadBuses(i)) * dp*Pd0(loadBuses(j));
                end
            end
        end
    end
    %for N_MC = 1000 : 1000 : 1e5
    for N_MC = 1e6

        iteste = iteste + 1;

        fprintf('Teste No %d - SMC \n', iteste);

        tic;
        
        run_MCS;

        toc;

        erroMC = 2.5758*std(p)/sqrt(N_MC)/mean(p);   %2.5758 = icdf('norm',(1-0.99)/2,0,1)

        resultadoFinal(iteste).N_MC = N_MC;
        resultadoFinal(iteste).MCS_7porc = p(end)*100;
        resultadoFinal(iteste).tempo = toc;
        resultadoFinal(iteste).erro = erroMC;
        resultadoFinal(iteste).txConverg = i/casos*100;        
     end
end
%save('resultadoFinal.mat','resultadoFinal','-mat');

%%
% xx = vertcat(resultadoFinal.N_MC);
% yy1 = vertcat(resultadoFinal.erro);
% yy2 = vertcat(resultadoFinal.tempo);
% yyaxis left
% plot(xx,yy1);
% grid minor;   
% xlabel('Número de amostras usados na SMC');
% ylabel('Erro');  
% grid minor;
% yyaxis right
% plot(xx,yy2);
% ylabel('Tempo');
% set(gca,'FontSize',20,'FontName','Times');
% set(gcf,'Paperunits','inches','PaperPosition',[0 0 20 10],'PaperSize',[20 10])
% pause(1);
% print('TempoSMC.pdf','-dpdf','-r400');

