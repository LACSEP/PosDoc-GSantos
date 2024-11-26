clc; clear all; 

addpath('.\Sistemas')  %Pasta com os sitemas testes

% casebase -> CARTAO NO FORMATO ADOTADO PELO MATHEUS!
casobase = case5;
[nbus,narea,area,pq,pv,ref,gen,Ybus,Vm,Va,Pg0,Qg,Pd0,Qd0,Qmax2,Qmin2,Vmax,Vmin,Sb,names,buses] = data(casobase);
Qmax = zeros(nbus,1); Qmax([ref;pv]) = Qmax2;
Qmin = zeros(nbus,1); Qmin([ref;pv]) = Qmin2;

%% ---- Ordering and size of data vectors ---- %%
npv = length(pv); npq = length(pq);         
gen = sort(gen); nger = length(gen);

Sd0 = sqrt(Pd0.^2+Qd0.^2); 
pos = find(Sd0==0); npos = length(pos);
pos2 = find(Sd0~=0); npos2 = length(pos2);
bp = Pd0./Sd0;  bp(ref) = 0; bp(pv) = 0;
bq = Qd0./Sd0;  bq(ref) = 0; bq(pv) = 0;

t = 0; %Load Growth Parameter
wP = ones(nbus,1);
wP(ref) = 0;
wQ = ones(nbus,1);
wG = ones(nger,1);
Vr = Vm;
Qspc = Qg(gen);
Qmax(ref) = 999.99;
Qmin(ref) = -99.99;
tol = 1e-4;   %Tolerance (Newton-Raphson)
itMax = 20;

%%%%%%%%%%%%%%%%%%%%Caso base%%%%%%%%%%%%%%%%%%%%%%
kpq = zeros(nbus,1); kpq(pq) = 1;
[tloadCB dPtotalCB it] = run_DirectMethod2(t, wP, wQ, tol, pq, npq, pv, npv,Ybus,Vm,Va,Pg0,Qg,Pd0,Qd0,nbus,itMax,kpq,bp,bq);
PdCasoBase = Pd0 + kpq*tloadCB.*bp;
fprintf('-----RESULTADO DO CASO BASE-----\n');
fprintf('Acrescimo de carga total igual a %f p.u. com lambda igual a %f \n', dPtotalCB, tloadCB);


%% ---- Execucao do algoritmo de entropia cruzada ---- %%
rng('shuffle'); %Garantir a aleatoriedade

%Parâmetros do algoritmo
nsamples_EC = 50; %Número de amostras
nElite = 10; %Número de amostras de elite
itMaxEC = 100; sig2_min = 1e-2; %Criterio de parada da entropia cruzada
t = 0.8; v = 0.2; %Fatores de suavização

%Parâmetros do sistema
mu = ones(nbus,1); %Partindo do cenário do caso base
sig2 = 0.1.*ones(nbus,1);
for i = 1 : 1 : npos
    mu(pos(i)) = 0;
    sig2(pos(i)) = 0;
end

%Exclusao de direcao de crescimento para as barras com carga
it_EC = 0; percViol = 4/100; atingiu = 0;

while true
    it_EC = it_EC + 1; %Atualizacao do numero de iteracoes
    %% Passo 1 - Criar o veotr de amostras X e avaliar a funcao objetivo fEC
    X = zeros(nbus,nsamples_EC);
    for i = 1 : 1 : nsamples_EC
        %Passo 1 - Criar o vetor de amostras X e avaliar a funcao objetivo
        while true
            X(:,i) = mu + sqrt(sig2.*abs(randn(nbus,1)));
            for j = 1 : 1 : npos
                X(pos(j),i) = 0;
            end
            kpq = X(:,1);
            [tload dPtotal it] = run_DirectMethod2(t, wP, wQ, tol, pq, npq, pv, npv, Ybus, Vm, Va, Pg0, Qg, Pd0, Qd0, nbus, itMax, kpq, bp, bq);
            fpen = 0;
            if tload > 0 && it <= itMax
                fEC(i,1) = dPtotal + fpen;
                fEC(i,2) = tload;
                fEC(i,3) = dPtotal;
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
        if  min(fEC(idxElite(i),2)) < percViol
            nviol = nviol + 1;
        end
    end
    if nviol >= 0.8*nElite && atingiu == 0
        atingiu = 1;
        fprintf('\n-----RESULTADO DO MÉTODO DE ENTROPIA-----\n');
        fprintf('Atingiu a condicao de violacao de tensao em %d iteracoes \n', it_EC);
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
    if it_EC == itMaxEC || max(sig2) < sig2_min
        if atingiu == 0
            fprintf('\n-----RESULTADO DO MÉTODO DE ENTROPIA-----\n');
        end
        [~,idx] = min(fEC(:,1));
        fprintf('Resultado final da entropia \n')
        fprintf('Acrescimo de carga total igual a %f p.u. com lambda igual a %f \n', fEC(idx,1), fEC(idx,2));
        fprintf('Valor de kp*lambda:\n');
        for i = 1 : 1 : nbus
            fprintf('Barra: %d: %.4f \n', buses(i),fEC(idx,2)*X(i,idx).*bp(i));
        end
        for i = 1 : 1 : nElite
            for j = 1 : 1 : nbus
                PdElite(i,j) = Pd0(j) + fEC(idxElite(i),2)*Xel(j,i)*bp(j);
            end
        end
        break;
    end
end
[~,idxElite] = sort(fEC(:,1));
Pd_EC = Pd0 + X(:,idx)*fEC(idx,2).*bp;

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

it_ED = 0;
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
                else
                    novapop(pos2(j),i) = pop(pos2(j),i);
                end
            end
        elseif estrategia(i) == 4
            F2 = rand(); %Fator de mutação 0 a 1
            for j = 1 : 1 : npos2
                novapop(j,i) = pop(pos2(j),i) + F2 * (pop(pos2(j),sel(1)) - pop(pos2(j),i)) + F * (pop(pos2(j),sel(2))-pop(pos2(j),sel(3)));
            end
        end
        kpq = novapop(:,i);
        [tload dPtotal it] = run_DirectMethod2(t, wP, wQ, tol, pq, npq, pv, npv, Ybus, Vm, Va, Pg0, Qg, Pd0, Qd0, nbus, itMax, kpq, bp, bq);
        fpen = 0;
        fnovapop(i,1) = dPtotal;
        fnovapop(i,2) = tload;
        fnovapop(i,3) = dPtotal;
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
        fprintf('Acrescimo de carga total igual a %f p.u. com lambda igual a %f \n', fpop(idx,1), fpop(idx,2));
        fprintf('Valor de kp*lambda:\n');
        for i = 1 : 1 : nbus
            fprintf('Barra: %d: %.4f \n', buses(i),fpop(idx,2)*pop(i,idx).*bp(i));
        end
        break;
    end
end
Pd_ED = Pd0+fpop(idx,2)*pop(:,idx).*bp;

% 
%Plotagem do grafico 
nc = 0;
for i = 0 : 0.1 : 1    
    for j = 0 : 0.1 : 1
        for k = 0 : 0.1 : 1 
            if i ~= 0 || j ~= 0 || k ~= 0
                kpq(pos2(1)) = i;
                kpq(pos2(2)) = j;
                kpq(pos2(3)) = k;
                [tload dPtotal it] = run_DirectMethod2(t, wP, wQ, tol, pq, npq, pv, npv,Ybus,Vm,Va,Pg0,Qg,Pd0,Qd0,nbus,itMax,kpq,bp,bq);
                if isnan(tload) == 0 && tload > 0 && it <= itMax 
                    nc = nc + 1;
                    res(nc).t = tload;
                    res(nc).dPload = dPtotal; 
                    res(nc).kpq2 = kpq(2); 
                    res(nc).kpq4 = kpq(4); 
                    res(nc).kpq5 = kpq(5); 
                    res(nc).Pd2 = Pd0(2) + kpq(2)*tload*bp(2); 
                    res(nc).Pd4 = Pd0(4) + kpq(4)*tload*bp(4); 
                    res(nc).Pd5 = Pd0(5) + kpq(5)*tload*bp(5);                    
                end
            end
        end        
    end
end

x = vertcat(res.Pd2);
y = vertcat(res.Pd4);
z = vertcat(res.Pd5);
npts = length(x);
% for i = 1 : 1 : npts
%     for j = 1 : 1 : npts
%         Z(i,j) = z(i);
%     end
% end
scatter3(x,y,z,'.');
hold on;
axis([0 max(x) 0 max(y) 0 max(z)]);
x2 = []; y2 = []; z2 = [];
xlabel('Pl_2'); ylabel('Pl_4'); zlabel('Pl_5');
for i = 1 : 1 : npts
    if x(i) == Pd0(2)
        scatter3(x(i),y(i),z(i),'b.');
        x2 = [x2 [x(i)]];
        y2 = [y2 [y(i)]];
        z2 = [z2 [z(i)]];
    elseif y(i) == Pd0(4)
        scatter3(x(i),y(i),z(i),'b.');
        x2 = [x2 [x(i)]];
        y2 = [y2 [y(i)]];
        z2 = [z2 [z(i)]];
    elseif z(i) == Pd0(5)
        scatter3(x(i),y(i),z(i),'b.');
        x2 = [x2 [x(i)]];
        y2 = [y2 [y(i)]];
        z2 = [z2 [z(i)]];
    else
    %    scatter3(x(i),y(i),z(i),'r.');
    end
end
hold on
%Caso base
scatter3(Pd0(2),Pd0(4),Pd0(5),'r');
scatter3(Pd0(2)+bp(2)*tloadCB,Pd0(4)+bp(4)*tloadCB,Pd0(5)+bp(5)*tloadCB,'k');
plot3([Pd0(2),Pd0(2)+bp(2)*tloadCB],[Pd0(4),Pd0(4)+bp(4)*tloadCB],[Pd0(5),Pd0(5)+bp(5)*tloadCB],'k');
grid minor;

for i = 1 : 1 : nElite
    plot3([Pd0(2),PdElite(i,2)],[Pd0(4),PdElite(i,4)],[Pd0(5),PdElite(i,5)],'r--');
    scatter3(PdElite(i,2),PdElite(i,4),PdElite(i,5),'r');
end

plot3([Pd0(2),Pd_ED(2)],[Pd0(4),Pd_ED(4)],[Pd0(5),Pd_ED(5)],'m--');
scatter3(Pd_ED(2),Pd_ED(4),Pd_ED(5),'m');

text(PdCasoBase(2),PdCasoBase(4),PdCasoBase(5),'\leftarrow On-scenario','FontSize',15,'FontName','Times')
text(Pd_EC(2),Pd_EC(4),Pd_EC(5),'\leftarrow Minimum-\lambda of cross-entropy algorithm','FontSize',15,'FontName','Times')
text(Pd_ED(2),Pd_ED(4),Pd_ED(5),'\leftarrow Minimum-\lambda of SADE algorithm','FontSize',15,'FontName','Times')
xlabel ('PL2(pu)','FontSize',15,'FontName','Times');
ylabel ('PL4(pu)','FontSize',15,'FontName','Times');
zlabel ('PL5(pu)','FontSize',15,'FontName','Times');


function [sel] = selecaoInd(nInd, targetInd)

aux = [1:nInd];
aux(targetInd) = [];

for i = 1 : 1 : 5
    k = randi ([1 nInd-i]);
    pos = aux(k);
    aux(k) = [];
    sel(i) = pos;
end

end
% 
% 
% set(gca,'FontSize',15,'FontName','Times')
% set(gcf,'Paperunits','inches','PaperPosition',[0 0 20 15],'PaperSize',[20 15])
% pause(1);
% print('Caso5','-dpdf','-r400')