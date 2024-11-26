clc; clear;

%Referência base: https://doi.org/10.1109/PESGM52003.2023.10252463

addpath('.\Sistemas')  %Pasta com os sitemas testes


%% ---- Caso teste---- %%
% casebase -> CARTAO NO FORMATO ADOTADO PELO MATHEUS!
casobase = case3_mod;
%casobase = case14;
%casobase = case107;
%casobase = case118
[nbus,narea,area,pq,pv,ref,gen,Ybus,Vm,Va,Pg0,Qg,Pd0,Qd0,Qmax,Qmin,Vmax,Vmin,Sb,names,buses] = data(casobase);
npv = length(pv); npq = length(pq);         %Vector Length

%% ---- Parametros iniciais do metodo direto ---- %%
t = 0;                                      %Load Growth Parameter
wP = zeros(size(Vm)); wQ = zeros(size(Vm));
wP([pv; pq]) = 1; wQ(pq) = 1;                %Initial Eigenvector
tol = 1e-4;   %Tolerance (Newton-Raphson)
itMax = 20;

%% ---- Execucao do metodo direto para o caso base ---- %%
[tload dPtotal it] = run_DirectMethod(t, wP, wQ, tol, pq, npq, pv, npv, Ybus, Vm, Va, Pg0, Qg, Pd0, Qd0, nbus, itMax);
PdCasoBase = (1+tload)*Pd0;
fprintf('-----RESULTADO DO CASO BASE-----\n');
fprintf('Acrescimo de carga total igual a %f p.u. com lambda igual a %f \n', dPtotal, tload);

%% ---- Execucao do algoritmo de entropia cruzada ---- %%
rng('shuffle'); %Garantir a aleatoriedade

%Parâmetros do algoritmo
nsamples_EC = 50; %Número de amostras
nElite = 10; %Número de amostras de elite
t = 0.8; v = 0.2; %Fatores de suavização
%itMaxEC = 100;
itMaxEC = 200;
%Media e desvio padrao das barras com carga
pos = find(Pd0 == 0); npos = length(pos);
pos2 = find(Pd0 ~= 0); npos2 = length(pos2);
[mu_load sig2_load] = loadData('case3mod');
a = 1./Pd0; b = -1;
%mu = a.*mu_load+b; mu(pos) = 0;
mu = [0;1;1]; %Partindo do cenário do caso base
muRef = mu;
sig2 = (a.^2).*sig2_load; sig2(pos) = 0;

it_EC = 0; atingiu = 0;
fprintf('\n-----RESULTADO DO MÉTODO DE ENTROPIA-----\n');
while true   
    it_EC = it_EC + 1; %Atualizacao do numero de iteracoes   
    %% Passo 1 - Criar o veotr de amostras X e avaliar a funcao objetivo f_EC
    X = zeros(nbus,nsamples_EC);
    for i = 1 : 1 : nsamples_EC
        while true
            X(:,i) = mu + sqrt(sig2.*abs(randn(nbus,1)));
            [SX dPtotal it] = run_DirectMethod(t, wP, wQ, tol, pq, npq, pv, npv, Ybus, Vm, Va, Pg0, Qg, Pd0, Qd0, nbus, itMax, X(:,i));
            fpen = 0;
%             for j = 1 : 1 : npos2
%                 fpen = fpen + (X(pos2(j))-muRef(pos2(j)))^2*(1+Pd0(pos2(j)));
%             end
           % fpen = min(10*sum(muRef),fpen);     
            if SX > 0 && it <= itMax 
                f_EC(i,1) = dPtotal + fpen;
                f_EC(i,2) = SX;
                f_EC(i,3) = dPtotal;
                break;
            end
        end
    end
    
    %% Passo 2 - Criar o conjunto ordenado por S(X) que contém as amostras de elite
    [~,idx_Elite] = sort(f_EC(:,1));
    Xel = X(:,idx_Elite(1:nElite)); 
  
    %% Passo 3 - Identificar se ha violacao de margem de seguranca
    nviol = 0;
    for i = 1 : 1 : nElite
        if  min(f_EC(idx_Elite(i),2)) < 0.01
            nviol = nviol + 1;
        end
    end
    if nviol >= 0.8*nElite && atingiu == 0
        atingiu = 1;
        fprintf('Atingiu a condicao de violacao de tensao em %d iteracoes \n', it_EC);
        fprintf('Carga 2: mu = %f e sigma = %f \n', mu(2), sig2(2));
        fprintf('Carga 3: mu = %f e sigma = %f \n', mu(3), sig2(3));
        muViolEC = mu;
        sig2ViolEC = sig2;
    end
    
    %% Passo 4 - Recalcular a média e a variância
    mu_0 = mu; sig2_0 = sig2;
    mu = mean(Xel')';
    sig2 = var(Xel')';
   
    %% Passo 5 - Aplicar fatores de suavização
    mu = t*mu+(1-t)*mu_0;
    sig2 = v*sig2+(1-v)*sig2_0;

    %% Passo 6 - Verificar condição de parada    
    if it_EC == itMaxEC || max(sig2) < 1e-3
        [~,idx] = min(f_EC(:,1));
        fprintf('Resultado final da entropia \n')
        fprintf('Acrescimo de carga total igual a %f p.u. com lambda igual a %f \n', f_EC(idx,1), f_EC(idx,2));
        fprintf('Valor de kp*lambda:\n');        
        for i = 1 : 1 : nbus
            fprintf('Barra: %d: %.4f \n', buses(i),f_EC(idx,2)*X(i,idx));
        end
        for i = 1 : 1 : nElite
            Pd_el(i,1) = (1 + f_EC(idx_Elite(i),2)*Xel(2,i))*Pd0(2);
            Pd_el(i,2) = (1 + f_EC(idx_Elite(i),2)*Xel(3,i))*Pd0(3);
        end
        break;
    end
end
[~,idx_Elite] = sort(f_EC(:,1));
Pd_EC = Pd0 + X(:,idx)*f_EC(idx,2).*Pd0;
  
% N1 = 1000; cont = 0;
% for i = 1 : 1 : N1    
%     kps = normrnd(muViolEC,sqrt(sig2ViolEC));   
%     f = normpdf(kps,mu_0,sqrt(sig2_0));
%     f(1) = 0;
%     q = normpdf(kps,muViolEC,sqrt(sig2ViolEC));
%     q(1) = 0;
%     [SX dPtotal it] = run_DirectMethod(t, wP, wQ, tol, pq, npq, pv, npv, Ybus, Vm, Va, Pg0, Qg, Pd0, Qd0, nbus, itMax, kps);
%     if dPtotal <  0.01
%         cont = cont + (f/g)
%     end
% end
% cont = cont/N1;


%% ---- Execucao do algoritmo de evolucao diferencial ---- %%
%Populacao do ED
pop = X(:,idx_Elite(1:nElite)); fpop = f_EC(idx_Elite(1:nElite),:);
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
    if it_ED == 0
        novapop = pop;
        fnovapop = fpop;
    end
    [melhorInd posMelhor] = min(fpop(:,1));
    it_ED = it_ED + 1;
    for i = 1 : 1 : nInd
        F = normrnd(0.5,0.1); %Gera um valor da gaussiana com média 0,5 e desvio padrão 0,1
        sel = selecaoInd(nInd, i);
        if it_ED <= LP
            estrategia(i) = randi(k);  %Durante o período de treinamento todas as estratégias têm a mesma chance de serem escolhidas
        else
            estrategia(i) = estrategia_final;
        end
        jrand = randi(pos2); %Obrigatoriamente um elemento de cada indivíduo sofrerá mutação
        if estrategia(i) == 1
            while true
                pc = normrnd(CRm(1), 0.1);
                if pc > 0 && pc <1
                    break;
                end
            end
            CRMemory(i) = pc;
            for j = 1 : 1 : npos2
                if rand() < pc || j == rand()
                    novapop(pos2(j),i) = pop(pos2(j),sel(1)) + F*pop(pos2(j),sel(2))-pop(pos2(j),sel(3));
                else
                    novapop(pos2(j),i) = pop(pos2(j),i);
                end
            end
        elseif estrategia(i) == 2
            while true
                pc = normrnd(CRm(1), 0.1);
                if pc > 0 && pc<1
                    break;
                end
            end
            CRMemory(i) = pc;
            for j = 1 : 1 : npos2
                if rand() < pc || j == rand()
                    novapop(pos2(j),i) = pop(pos2(j),i) + F*(pop(pos2(j),posMelhor)-pop(pos2(j),i)) + F*(pop(pos2(j),sel(1))-pop(pos2(j),sel(2))) + F*(pop(pos2(j),sel(3))-pop(pos2(j),sel(4)));
                else
                    novapop(pos2(j),i) = pop(pos2(j),i);
                end
            end
        elseif estrategia(i) == 3
            while true
                pc = normrnd(CRm(1), 0.1);
                if pc > 0 && pc <1
                    break;
                end
            end
            CRMemory(i) = pc;
            for j = 1 : 1 : npos2
                if rand() < pc || j == rand()
                    novapop(pos2(j),i) = pop(pos2(j),i) + F*pop(pos2(j),sel(2))-pop(pos2(j),sel(3))  + F*pop(pos2(j),sel(4))-pop(pos2(j),sel(5));
                else
                    novapop(pos2(j),i) = pop(pos2(j),i);
                end
            end
        elseif estrategia(i) == 4
            F2 = rand(); %Fator de mutação 0 a 1
            for j = 1 : 1 : npos2
                novapop(j,i) = novapop(pos2(j),i) + F2 * (pop(pos2(j),sel(1)) - pop(pos2(j),i)) + F * (pop(pos2(j),sel(2))-pop(pos2(j),sel(3)));
            end
        end
        [SX dPtotal it] = run_DirectMethod(t, wP, wQ, tol, pq, npq, pv, npv, Ybus, Vm, Va, Pg0, Qg, Pd0, Qd0, nbus, itMax, novapop(:,i));
        fpen = 0;
        for j = 1 : 1 : npos2
            fpen = fpen + (novapop(pos2(j))-muRef(pos2(j)))^2*(1+Pd0(pos2(j)));
        end
        %fpen = max(10*sum(muRef),fpen);
        fnovapop(i,1) = dPtotal;
        fnovapop(i,2) = SX;
        fnovapop(i,3) = dPtotal;
        %         j = 0;
        %         while true
        %             j = j + 1;
        %             if abs((mu(pos2(j),1) - muRef(pos2(j),1))/muRef(pos2(j),1)) > 0.5
        %                 crcruzamento = 1;
        %                 break;
        %             end
        %             if j == npos2
        %                 break;
        %             end
        %         end
        %   if crcruzamento == 1 || SX < 0 || it >= itMax || min(novapop(:,i)) < 0
        if  SX < 0 || it >= itMax || min(novapop(:,i)) < 0
            novapop(:,i) = pop(:,i);
            fnovapop = fpop;
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
            taxa_sucesso(i) = ns(i) / (ns(i)+nf(i)) + cte;
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
            fprintf('Barra: %d: %.4f \n', buses(i),fpop(idx,2)*pop(i,idx));
        end
        break;
    end
end
Pd_ED = (1+fpop(idx,2)*pop(:,idx)).*Pd0;

%Plotagem dos graficos
% Kp(1,1) = 0;
% k = 0;
% for i = 0 : 0.01 : 1.2
%     Kp(2,1) = i;
%     for j = 0 : 0.01 : 1.2
%         Kp(3,1) = j;
%         if sum(Kp) ~= 0
%             [tload,dPtotal,it] = run_DirectMethod(t, wP, wQ, tol, pq, npq, pv, npv, Ybus, Vm, Va, Pg0, Qg, Pd0, Qd0, nbus, itMax, Kp);
%             if isnan(tload) == 0 && tload > 0 && it <= itMax
%                 k = k + 1;
%                 res(k).t = tload;
%                 res(k).dPtotal = dPtotal;
%                 res(k).Kp2 = Kp(2);
%                 res(k).Kp3 = Kp(3);
%                 res(k).Pd2 = (1+Kp(2)*tload)*Pd0(2);
%                 res(k).Pd3 = (1+Kp(3)*tload)*Pd0(3);
%             end
%         end
%     end
% end
% x = vertcat(res.Pd2);
% y = vertcat(res.Pd3);
x = load('res_x2.mat');
x = x.x;
y = load('res_y2.mat');
y = y.y;
scatter(x,y);
%Caso base
hold on;
scatter(Pd0(2),Pd0(3),'b')
plot([Pd0(2),PdCasoBase(2)],[Pd0(3),PdCasoBase(3)],'k--');
scatter(PdCasoBase(2),PdCasoBase(3),'k')
plot([Pd0(2),Pd_EC(2)],[Pd0(3),Pd_EC(3)],'r--');
scatter(Pd_EC(2),Pd_EC(3),'b')
plot([Pd0(2),Pd_ED(2)],[Pd0(3),Pd_ED(3)],'m--');
scatter(Pd_ED(2),Pd_ED(3),'m')
% 
% xlim([0 1.2])
% ylim([0 2.2])

for i = 1 : 1 : nElite
    plot([Pd0(2),Pd_el(i,1)],[Pd0(3),Pd_el(i,2)],'r--');
    scatter(Pd_el(i,1),Pd_el(i,2),'r');
end

text(PdCasoBase(2),PdCasoBase(3),'\leftarrow On-scenario','FontSize',15,'FontName','Times')
text(Pd_EC(2),Pd_EC(3),'\leftarrow Minimum-\lambda of cross-entropy algorithm','FontSize',15,'FontName','Times')
text(Pd_ED(2),Pd_ED(3),'\leftarrow Minimum-\lambda of SADE algorithm','FontSize',15,'FontName','Times')
xlabel ('P_2(pu)','FontSize',15,'FontName','Times');
ylabel ('P_3(pu)','FontSize',15,'FontName','Times');

grid minor;

% set(gca,'FontSize',15,'FontName','Times')
% set(gcf,'Paperunits','inches','PaperPosition',[0 0 20 15],'PaperSize',[20 15])
% pause(1);
% print('Caso2','-dpdf','-r400')




