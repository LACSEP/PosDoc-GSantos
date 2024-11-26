clc; clear all;

tic;

addpath('.\Sistemas')  %Pasta com os sitemas testes

% casebase -> CARTAO NO FORMATO ADOTADO PELO MATHEUS!
casobase = case5;
[nbus,narea,area,pq,pv,ref,gen,Ybus,Vm,Va,Pg0,Qg,Pd0,Qd0,Qmax2,Qmin2,Vmax,Vmin,Sb,names,buses] = data(casobase);
Qmax = zeros(nbus,1); Qmax(gen) = Qmax2;
Qmin = zeros(nbus,1); Qmin(gen) = Qmin2;

%% ---- Ordering and size of data vectors ---- %%
npv = length(pv); npq = length(pq);
gen = sort(gen); nger = length(gen);

%% ---- Initial Values ---- %%
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

%% ---- Base case power flow ---- %%
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
    f = [eP(1:end ~= ref);eQ;ey];  %Function f
    s = mldivide(df_dx,f);      %Convergence Error
    x = [Va(1:end ~= ref); Vm; Qspc];
    x = x-s;      %Next State Vector
    Va(1:end ~= ref) = x(1:nbus-1);
    Vm = x(nbus:2*nbus-1);
    Qspc = x(2*nbus:2*nbus+npv);
end

Sd0 = sqrt(Pd0.^2+Qd0.^2);
Pd0_total = sum(Pd0);
pos = find(Sd0==0); npos = length(pos);
pos2 = find(Sd0~=0); npos2 = length(pos2);
alphaq = Qd0./Pd0;  alphaq(pos) = 0;
Pg0(gen) = (real(Scalc(gen))-Pd0(gen));
fpart = zeros(nbus,1);
fpart(gen) = Pg0(gen)/sum(Pg0);
%% ---- Initial Values ---- %%
t = 0; %Load Growth Parameter
wP = ones(nbus,1);
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
tol = 1e-7;   %Tolerance (Newton-Raphson)
itMax = 200;

%%%%%%%%%%%%%%%%%%%%Caso base%%%%%%%%%%%%%%%%%%%%%%
bp = zeros(nbus,1); bp = Pd0;
[tloadCB,dPtotalCB,it] = run_DirectMethodQlim(t, wP, wQ, wG, tol, gen, nger, pq, npq, ref, pv, npv,Ybus,Vm,Va,Pg0,Qg,Qmin,Qmax,Pd0,nbus,itMax,bp,alphaq,fpart);
PdCasoBase = Pd0 + tloadCB.*bp;
fprintf('-----RESULTADO DO CASO BASE-----\n');
fprintf('Acrescimo de carga total igual a %f p.u. (%.2f%%) \n', dPtotalCB, dPtotalCB/Pd0_total*100);


%% ---- Execucao do algoritmo de entropia cruzada ---- %%
rng('shuffle'); %Garantir a aleatoriedade

%Parâmetros do algoritmo
nsamples_EC = 100; %Número de amostras
nElite = 90; %Número de amostras de elite
itMaxEC = 1000; sig2_min = 1e-4; %Criterio de parada da entropia cruzada
t = 0.8; v = 0.2; %Fatores de suavização

%Parâmetros do sistema
for i = 1 : 1 : npos2
    mu(pos2(i),1) = Pd0(pos2(i));  %Partindo do cenário do caso base
    sig2(pos2(i),1) = sqrt(0.1*Pd0(pos2(i)));
end

%Exclusao de direcao de crescimento para as barras com carga
it_EC = 0; percViol = 4/100; atingiu = 0;

while true
    it_EC = it_EC + 1; %Atualizacao do numero de itera coes
    %% Passo 1 - Criar o veotr de amostras X e avaliar a funcao objetivo fEC
    X = zeros(nbus,nsamples_EC);
    for i = 1 : 1 : nsamples_EC
        %Passo 1 - Criar o vetor de amostras X e avaliar a funcao objetivo
        while true
            for j = 1 : 1 : npos2
                while true 
                    X(pos2(j),i) = normrnd(mu(pos2(j)),sqrt(sig2(pos2(j))));
                    if X(pos2(j),i) > 0 
                        break;
                    end
                end
            end
            bp = X(:,i);
            [tload,dPtotal,it] = run_DirectMethodQlim(t, wP, wQ, wG, tol, gen, nger, pq, npq, ref, pv, npv,Ybus,Vm,Va,Pg0,Qg,Qmin,Qmax,Pd0,nbus,itMax,bp,alphaq,fpart);
            fpen = 0;
            if tload > 0 && it <= itMax
                fEC(i,1) = dPtotal + fpen;
                fEC(i,2) = tload;
                fEC(i,3) = dPtotal;
                fEC(i,4) = dPtotal/Pd0_total; %Resultado para o artigo!
                break;
            end
        end
    end

    %% Passo 2 - Criar o conjunto ordenado por S(X) que contém as amostras de elite
    [~,idxElite] = sort(fEC(:,1));
    Xel = X(:,idxElite(1:nElite));

    %% Passo 3 - Identificar se ha violacao de margem de seguranca
    nviol = 0;
    for i = 1 : 1 : nElite
        if  fEC(idxElite(i),4) < percViol
            nviol = nviol + 1;
        end
    end
    if nviol >= 0.8*nElite && atingiu == 0
        atingiu = 1;
        fprintf('\n-----RESULTADO DO MÉTODO DE ENTROPIA-----\n');
        fprintf('Atingiu a condicao de violacao de margem de estabilidade em %d iteracoes \n', it_EC);
        muViolEC = mu;
        sig2ViolEC = sig2;
    end

    %% Passo 4 - Recalcular a média e a variância
    mu_0 = mu; sig2_0 = sig2;
    for i = 1 : 1 : nbus
        mu(i,1) = mean(Xel(i,:));
        sig2(i,1) = var(Xel(i,:));
    end

    %% Passo 5 - Aplicar fatores de suavização
    mu = t*mu+(1-t)*mu_0;
    sig2 = v*sig2+(1-v)*sig2_0;

    %% Passo 6 - Verificar condição de parada
    if it_EC == itMaxEC %|| max(sig2) < sig2_min
        if atingiu == 0
            fprintf('\n-----RESULTADO DO MÉTODO DE ENTROPIA-----\n');
        end
        [~,idx] = min(fEC(:,1));
        fprintf('Resultado final da entropia \n')
        fprintf('Acrescimo de carga total igual a %f p.u. (%.2f%%) \n', fEC(idx,1), fEC(idx,4)*100);
        fprintf('Valor de bp*lambda:\n');
        for i = 1 : 1 : nbus
            fprintf('Barra: %d: %.4f \n', buses(i),fEC(idx,2)*X(i,idx));
        end
        for i = 1 : 1 : nElite
            for j = 1 : 1 : nbus
                PdElite(i,j) = Pd0(j) + fEC(idxElite(i),2)*Xel(j,i);
            end
        end
        break;
    end
end
[~,idxElite] = sort(fEC(:,1));
Pd_EC = Pd0 + X(:,idx)*fEC(idx,2);
toc;
%% Passo 7 - Probabilidade de violacao da margem
if atingiu == 1
    N = 10000;
    c = 0;
    %Parâmetros do sistema
    for i = 1 : 1 : npos2
        mu(pos2(i),1) = Pd0(pos2(i));  %Partindo do cenário do caso base
        sig2(pos2(i),1) = sqrt(0.1*Pd0(pos2(i)));
    end
    i = 0;
    while true
        i = i + 1;
        while true
            bp = normrnd(muViolEC,sqrt(sig2ViolEC));
            [tload,dPtotal,it] = run_DirectMethodQlim(t, wP, wQ, wG, tol, gen, nger, pq, npq, ref, pv, npv,Ybus,Vm,Va,Pg0,Qg,Qmin,Qmax,Pd0,nbus,itMax,bp,alphaq,fpart);
            if tload > 0 && it <= itMax
                fMC(i,1) = dPtotal/Pd0_total;
                break;
            end
        end
        if fMC(i) < percViol
            p = 1; g = 1;
            for j = 1 : 1 : npos2
                pf(j) = 1/sqrt(2*pi*sig2(pos2(j)))*exp(-(bp(pos2(j))-mu(pos2(j)))^2/(2*sig2(pos2(j))));
                pg(j) = 1/sqrt(2*pi*sig2ViolEC(pos2(j)))*exp(-(bp(pos2(j))-muViolEC(pos2(j)))^2/(2*sig2ViolEC(pos2(j))));
                p = p*pf(j);
                g = g*pg(j);
            end
            c = c + p/g;
        end
        if i > N && abs(fMC(i) - fMC(i-500)) < 10^-6
            fprintf('Probabilidade de violacao para %d iteracoes: %.4f%% \n', i, c/i*100);
            break;
        end
    end
end


%% ---- Execucao do algoritmo de evolucao diferencial ---- %%
%Populacao do ED
pop = X(:,idxElite(1:nElite)); fpop = fEC(idxElite(1:nElite),:);
nInd = nElite;

itMaxED = 100;
LP = 20; %Periodo de aprendizagem
k = 4; %Quantidade de mutacoes avaliadas
cont = zeros(1,k);
CRm = ones(1,k)*0.5;  %Inicializa a taxa de crossover
p = 1/k*ones(1,k); %Inicializa as probabilidades de escolher cada mutacao
ns = zeros(1,k); %Numero de mutacoes bem sucedidas por cada uma das estrategias
nf = zeros(1,k); %Numero de falhas de cada uma das estrategias

it_ED = 0; t = 0;
while true
    it_ED = it_ED + 1;
    if it_ED == 1
        novapop = pop;
        fnovapop = fpop;
    end
    [melhorInd posMelhor] = min(fpop(:,1));
    for i = 1 : 1 : nInd
        F = normrnd(0.5,0.1); %Gera um valor da gaussiana com média 0,5 e desvio padrão 0,1
        sel = selecaoInd(nInd, i);
        if it_ED <= LP
            estrategia(i) = randi(k);  %Durante o período de treinamento todas as estratégias têm a mesma chance de serem escolhidas
        else
            estrategia(i) = estrategia_final;
        end
        jrand = randi(npos2); %Obrigatoriamente um elemento de cada indivíduo sofrerá mutação
        if estrategia(i) == 1
            while true
                pc = normrnd(CRm(1), 0.1);
                if pc > 0 && pc <1
                    break;
                end
            end
            CRMemory(i) = pc;
            for j = 1 : 1 : npos2
                if rand() < pc || j == jrand
                    novapop(pos2(j),i) = pop(pos2(j),sel(1)) + F*(pop(pos2(j),sel(2))-pop(pos2(j),sel(3)));
                    if novapop(pos2(j),i) < 0
                        novapop(pos2(j),i) = pop(pos2(j),i);
                    end
                else
                    novapop(pos2(j),i) = pop(pos2(j),i);
                end
            end
        elseif estrategia(i) == 2
            while true
                pc = normrnd(CRm(2), 0.1);
                if pc > 0 && pc<1
                    break;
                end
            end
            CRMemory(i) = pc;
            for j = 1 : 1 : npos2
                if rand() < pc || j == jrand
                    novapop(pos2(j),i) = pop(pos2(j),i) + F*(pop(pos2(j),posMelhor)-pop(pos2(j),i)) + F*(pop(pos2(j),sel(1))-pop(pos2(j),sel(2))) + F*(pop(pos2(j),sel(3))-pop(pos2(j),sel(4)));
                    if novapop(pos2(j),i) < 0
                        novapop(pos2(j),i) = pop(pos2(j),i);
                    end
                else
                    novapop(pos2(j),i) = pop(pos2(j),i);
                end
            end
        elseif estrategia(i) == 3
            while true
                pc = normrnd(CRm(3), 0.1);
                if pc > 0 && pc <1
                    break;
                end
            end
            CRMemory(i) = pc;
            for j = 1 : 1 : npos2
                if rand() < pc || j == jrand
                    novapop(pos2(j),i) = pop(pos2(j),i) + F*pop(pos2(j),sel(2))-pop(pos2(j),sel(3))  + F*pop(pos2(j),sel(4))-pop(pos2(j),sel(5));
                    if novapop(pos2(j),i) < 0
                        novapop(pos2(j),i) = pop(pos2(j),i);
                    end
                else
                    novapop(pos2(j),i) = pop(pos2(j),i);
                end
            end
        elseif estrategia(i) == 4
            F2 = rand(); %Fator de mutação 0 a 1
            for j = 1 : 1 : npos2
                novapop(pos2(j),i) = pop(pos2(j),i) + F2 * (pop(pos2(j),sel(1)) - pop(pos2(j),i)) + F * (pop(pos2(j),sel(2))-pop(pos2(j),sel(3)));
                if novapop(j,i) < 0
                    novapop(j,i) = pop(pos2(j),i);
                end
            end
        end
        bp = novapop(:,i);
        [tload dPtotal it] = run_DirectMethodQlim(t, wP, wQ, wG, tol, gen, nger, pq, npq, ref, pv, npv,Ybus,Vm,Va,Pg0,Qg,Qmin,Qmax,Pd0,nbus,itMax,bp,alphaq,fpart);
        fpen = 0;
        fnovapop(i,1) = dPtotal;
        fnovapop(i,2) = tload;
        fnovapop(i,3) = dPtotal;
        fnovapop(i,4) = dPtotal/Pd0_total;
        %         if fnovapop(i,4) < 0.04
        %              pause()
        %         end
        if  tload < 0 || it >= itMax || min(novapop(:,i)) < 0
            novapop(:,i) = pop(:,i);
            fnovapop(i,:) = fpop(i,:);
        end
    end
    %%Processo de selecao
    for i = 1 : 1 : nInd
        if (fnovapop(i,1)) < fpop(i,1)
            pop(:,i) = novapop(:,i);
            fpop(i,:)= fnovapop(i,:);
            if it_ED <= LP
                ns(estrategia(i)) = ns(estrategia(i)) + 1;
                cont(estrategia(i)) = cont(estrategia(i)) + 1;
                if estrategia(i) ~= 4
                    CR(cont(estrategia(i)),estrategia(i)) = CRMemory(i);
                end
            end
        else
            if it_ED <= LP
                nf(estrategia(i)) = nf(estrategia(i)) + 1;
            end
        end
    end
    if it_ED == LP
        cte = 0.01;
        sumts = 0; %somatorio das taxas de sucesso
        estrategia_final = 1;
        for i = 1 : 1 : k
            taxa_sucesso(i) = ns(i) / (ns(i)+nf(i));
            sumts = sumts + taxa_sucesso(i);
        end
        for i = 1 : 1 : k
            prob(i) = taxa_sucesso(i) / sumts;
            if prob(i) > prob(estrategia_final)
                estrategia_final = i;
            end
            if i ~= 4
                CRm(i) = median(CR([1:cont(i)],i));
            end
        end
    end
    if it_ED == itMaxED
        [~,idx] = min(fpop(:,1));
        fprintf('\n-----RESULTADO DO MÉTODO DE EVOLUCAO DIFERENCIAL-----\n');
        fprintf('Acrescimo de carga total igual a %f p.u. (%.2f%%) \n', fpop(idx,1), fpop(idx,4)*100);
        fprintf('Valor de bp*lambda:\n');
        for i = 1 : 1 : nbus
            fprintf('Barra: %d: %.4f \n', buses(i),fpop(idx,2)*pop(i,idx));
        end
        break;
    end
end
Pd_ED = Pd0+fpop(idx,2)*pop(:,idx);

% %
% %Plotagem do grafico
% nc = 0; kpq = zeros(nbus,1);
% for i = 0 : 0.1 : 1
%     for j = 0 : 0.1 : 1
%         for k = 0 : 0.1 : 1
%             if i ~= 0 || j ~= 0 || k ~= 0
%                 kpq(pos2(1)) = i;
%                 kpq(pos2(2)) = j;
%                 kpq(pos2(3)) = k;
%                 [tload dPtotal it] = run_DirectMethodQlim(t, wP, wQ, wG, tol, gen, nger, pq, npq, ref, pv, npv,Ybus,Vm,Va,Pg0,Qg,Qmin,Qmax,Pd0,Qd0,nbus,itMax,kpq,bp,bq);
%                 if isnan(tload) == 0 && tload > 0 && it <= itMax
%                     nc = nc + 1;
%                     res(nc).t = tload;
%                     res(nc).dPload = dPtotal;
%                     res(nc).kpq2 = kpq(2);
%                     res(nc).kpq4 = kpq(4);
%                     res(nc).kpq5 = kpq(5);
%                     res(nc).Pd2 = Pd0(2) + kpq(2)*tload*bp(2);
%                     res(nc).Pd4 = Pd0(4) + kpq(4)*tload*bp(4);
%                     res(nc).Pd5 = Pd0(5) + kpq(5)*tload*bp(5);
%                 end
%             end
%         end
%     end
% end

% x = vertcat(res.Pd2);
% y = vertcat(res.Pd4);
% z = vertcat(res.Pd5);
% npts = length(x);
% % for i = 1 : 1 : npts
% %     for j = 1 : 1 : npts
% %         Z(i,j) = z(i);
% %     end
% % end
% scatter3(x,y,z,'.');
% hold on;
% axis([0 max(x) 0 max(y) 0 max(z)]);
% x2 = []; y2 = []; z2 = [];
% xlabel('Pl_2'); ylabel('Pl_4'); zlabel('Pl_5');
% for i = 1 : 1 : npts
%     if x(i) == Pd0(2)
%         scatter3(x(i),y(i),z(i),'b.');
%         x2 = [x2 [x(i)]];
%         y2 = [y2 [y(i)]];
%         z2 = [z2 [z(i)]];
%     elseif y(i) == Pd0(4)
%         scatter3(x(i),y(i),z(i),'b.');
%         x2 = [x2 [x(i)]];
%         y2 = [y2 [y(i)]];
%         z2 = [z2 [z(i)]];
%     elseif z(i) == Pd0(5)
%         scatter3(x(i),y(i),z(i),'b.');
%         x2 = [x2 [x(i)]];
%         y2 = [y2 [y(i)]];
%         z2 = [z2 [z(i)]];
%     else
%     %    scatter3(x(i),y(i),z(i),'r.');
%     end
% end
% hold on
% %Caso base
% scatter3(Pd0(2),Pd0(4),Pd0(5),'r');
% scatter3(Pd0(2)+bp(2)*tloadCB,Pd0(4)+bp(4)*tloadCB,Pd0(5)+bp(5)*tloadCB,'k');
% plot3([Pd0(2),Pd0(2)+bp(2)*tloadCB],[Pd0(4),Pd0(4)+bp(4)*tloadCB],[Pd0(5),Pd0(5)+bp(5)*tloadCB],'k');
% grid minor;
%
% for i = 1 : 1 : nElite
%     plot3([Pd0(2),PdElite(i,2)],[Pd0(4),PdElite(i,4)],[Pd0(5),PdElite(i,5)],'r--');
%     scatter3(PdElite(i,2),PdElite(i,4),PdElite(i,5),'r');
% end
%
% plot3([Pd0(2),Pd_ED(2)],[Pd0(4),Pd_ED(4)],[Pd0(5),Pd_ED(5)],'m--');
% scatter3(Pd_ED(2),Pd_ED(4),Pd_ED(5),'m');
%
% text(PdCasoBase(2),PdCasoBase(4),PdCasoBase(5),'\leftarrow On-scenario','FontSize',15,'FontName','Times')
% text(Pd_EC(2),Pd_EC(4),Pd_EC(5),'\leftarrow Minimum-\lambda of cross-entropy algorithm','FontSize',15,'FontName','Times')
% text(Pd_ED(2),Pd_ED(4),Pd_ED(5),'\leftarrow Minimum-\lambda of SADE algorithm','FontSize',15,'FontName','Times')
% xlabel ('PL2(pu)','FontSize',15,'FontName','Times');
% ylabel ('PL4(pu)','FontSize',15,'FontName','Times');
% zlabel ('PL5(pu)','FontSize',15,'FontName','Times');
%
%
% %%%%%Atencao, confirme o cartao do sistema!
% %save('res5busQlim','Pd0','PdCasoBase','PdElite','Pd_EC','Pd_ED','res')
% save('res5busQlimD','Pd0','PdCasoBase','PdElite','Pd_EC','Pd_ED','res')
%
% function [sel] = selecaoInd(nInd, targetInd)
%
% aux = [1:nInd];
% aux(targetInd) = [];
%
% for i = 1 : 1 : 5
%     k = randi ([1 nInd-i]);
%     pos = aux(k);
%     aux(k) = [];
%     sel(i) = pos;
% end
%
% end
% %
% %
% % set(gca,'FontSize',15,'FontName','Times')
% % set(gcf,'Paperunits','inches','PaperPosition',[0 0 20 15],'PaperSize',[20 15])
% % pause(1);
% % print('Caso5','-dpdf','-r400')