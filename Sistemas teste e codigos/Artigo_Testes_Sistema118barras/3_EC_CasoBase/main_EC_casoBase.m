clc; clear; 

warning off;

addpath('.\Sistemas')

%Sistema
system = 'case118bus_modificado.cdf';

%% ---- Leitura do sistema teste e obtenção da matriz de correlação ---- %%
run_initialization;
Pd0total = sum(Pd0(loadBuses));

%% ---- Método da entropia cruzada ---- %%

fprintf('----------------------MÉTODO DA EC----------------------\n')
%% Teste preliminar - obtenção da g(x)

%Parâmetros da entropia 
N_EC = 1000; %no. de amostras 
N_el = ceil(0.1*N_EC); %no. de elite 
itMax_EC = 100;

%Alvo 
percViol = 0.07;

tic; rng(1);

it_EC = 0;

while true 

    it_EC = it_EC + 1;    
    
    casos = 0;
    % --- Amostragem - Na primeira iteração, o sorteio é feito em X e nas demais em Z
    for i = 1 : 1 : N_EC        

        while true %Os sorteios sao repetidos ate o criterio de avaliacao ser satisfeito

            %Sorteio
            
            while true
                
                casos = casos + 1;

                if it_EC == 1 %Sorteio em X
                     X(i,:) = mvnrnd(mu,Sigma,1);
                
                else %Sorteio em Z
                    Z(i,:) = mvnrnd(muZ,sigmaZ,1);
                    for j = 1 : 1 : nLW
                        U(i,j) = normcdf(Z(i,j));
                    end
                    for j = 1 : 1 : nLW
                        Z(i,j) = norminv(U(i,j),0,1);
                    end
                    for j = 1 : 1 : nLW
                        X(i,j) = Z(i,j) * sig(j) + mu(j);
                    end
                end

                if min(X(i,100:106)) >= 0 %Garante que carga continua consumindo e eólica gerando
                    break;
                end

            end

            % Avaliação

            %1) Atualização dos valores de carga e geração
            bpl(loadBuses) = X(i,1:nLoad);
            Pd0(wind) = -1*X(i,nW1:end);
            Pd0total = sum(bpl);
            Pg0(1:end ~= ref) = Pg0_old(1:end ~= ref) + Pg0_old(1:end ~= ref)/Pgt_old*(2*7+sum(Pd0(wind)));

            %2) Cálculo da margem de carrregamento se o FP convergir            
            %Fluxo de carga
            [Vm, Va, Pg0, Qg, itPF] = runPowerFlow(Vm,Va,Vr,Pd0,Qd0,Pg0,Qg,nbus,npv,ref,gen,nger,Qg(gen),Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,Ybus,tepa,tepr,itPFmax);
            if itPF < itPFmax
                %Método Direto
                [tload, dPtotal, itPoC] = runDirectMethod(Vm,Va,Vr,Pd0,Pg0,Qg,nbus,npv,ref,gen,nger,Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,dA_dQg,dC_dVa,Ybus,bpl,tan0,fpart,tol,itPoCmax);
                %Verifica convergência do método direto
                if  itPoC <= itPoCmax
                    Pd = Pd0 + tload*bpl;
                    %Verifica se as cargas continuam maiores do que 0
                    if min(Pd(loadBuses)) >= 0
                        Pd = Pd0 + tload*bpl;
                        S_EC(i,1) = dPtotal/Pd0total;
                        S_EC(i,2) = tload;
                        S_EC(i,3) = Pd0total;
                        break;
                    end
                end
            end

        end
    end

    if it_EC == 1 

        %Normalize the original variables (mu = 0, std = 1)
        X_normalized = (X - mean(X)) ./ std(X);

        %Compute the empirical correlation matrix of the normalized variables
        empirical_R = corr(X);

        % Compute the Cholesky decomposition of the correlation matrix
        L_empirical = chol(empirical_R, 'lower');

        %
        L_target = chol(corrSigma, 'lower');
        
        % Adjust the normalized variables to match the target correlation structure
        Z = X_normalized / L_empirical' * L_target';

    end
    
    % --- Selecao das amostras de elite e atualização dos parametros

    %Sort in ascending order.
    [S,ix] = sort(S_EC(:,1));

    tr = S(N_el);

    if  tr < percViol || it_EC == itMax_EC
        fEC(:,1) = S_EC(ix,1);
        fEC(:,2) = S_EC(ix,2);
        fEC(:,3) = S_EC(ix,3);
        break;
    end

    if it_EC == 1
        fprintf('Resultado do método de EC: \n');
    end 
    
    fprintf('Iteracao %d da EC com %.2f%% dos casos convergidos: tr = %.4f \n', it_EC, i/casos*100, tr)

    %Elite sample
    elite_z = Z(ix(1:N_el),:);
    elite_x = X(ix(1:N_el),:);
    %Z-parameter update
    muZ = mean(elite_z);
    sigZ = std(elite_z);
    sigmaZ = cov(Z);
    %X-parameter update
    muX = mean(elite_x);
    sigX = std(elite_x);
    sigmaX = cov(X);   

end

%% Step 2 - MCS from CE results

N_EC2 = 1e4;

it_EC2 = 0;

I = []; W = []; p = []; 

rng(1); 

d = 10; casos = 0;

for i = 1 : 1 : N_EC2

    while true

        %Sorteio
        while true
            casos = casos + 1;
            X2(i,:) = mvnrnd(muX,sigmaX,1);
            if min(X2(i,100:106)) >= 0 %Garante que carga continua consumindo e eólica gerando
                break;
            end
        end

        % Avaliação

        %1) Atualização dos valores de carga e geração
        bpl(loadBuses) = X2(i,1:nLoad);
        Pd0(wind) = -1*X2(i,nW1:end);        
        Pg0(1:end ~= ref) = Pg0_old(1:end ~= ref) + Pg0_old(1:end ~= ref)/Pgt_old*(2*7+sum(Pd0(wind)));

        %2) Cálculo da margem de carrregamento se o FP convergir
        %Fluxo de carga
        [Vm, Va, Pg0, Qg, itPF] = runPowerFlow(Vm,Va,Vr,Pd0,Qd0,Pg0,Qg,nbus,npv,ref,gen,nger,Qg(gen),Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,Ybus,tepa,tepr,itPFmax);
        if itPF < itPFmax
            %Método Direto
            [tload, dPtotal, itPoC] = runDirectMethod(Vm,Va,Vr,Pd0,Pg0,Qg,nbus,npv,ref,gen,nger,Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,dA_dQg,dC_dVa,Ybus,bpl,tan0,fpart,tol,itPoCmax);
            %Verifica convergência do método direto
            if  itPoC <= itPoCmax
                Pd = Pd0 + tload*bpl;
                %Verifica se as cargas continuam maiores do que 0
                if min(Pd(loadBuses)) >= 0
                    Pd = Pd0 + tload*bpl;
                    S2_EC(i,1) = dPtotal/Pd0total;
                    S2_EC(i,2) = tload;
                    S2_EC(i,3) = Pd0total;
                    break;
                end
            end
        end
    end

    if S2_EC(i,1) < percViol
        I(i) = 1;
        f(i,1) = mvnpdf(X2(i,:),mu,Sigma);
        g(i,1) = mvnpdf(X2(i,:),muX,sigmaX);
        W(i) = f(i,1)/g(i,1);
    else
        I(i) = 0;
        W(i) = 0;
    end
    
    if sum(I) == 0
        p(i) = 0;
    else
        p(i) = 1/i*sum(I.*W);
    end

    if rem(i,d) == 0
        if d == 10
            fprintf('----------------------Resultados parciais do método de EC----------------------\n')
        end
        fprintf('Para %d iterações: %.4f%% com %.2f dos casos convergidos \n', d, p(end)*100, i/casos*100);
        d = d * 10;
    end
      
end

erroEC = 2.5758*std(p)/sqrt(N_EC2)/mean(p);   %2.5758 = icdf('norm',(1-0.99)/2,0,1)

toc;

resultadoFinal(1).EC_preAmost = N_EC;
resultadoFinal(1).EC_amost = N_EC2;
resultadoFinal(1).EC_4porc = p(end)*100;
resultadoFinal(1).tempo = toc;
resultadoFinal(1).erro = erroEC;
resultadoFinal(1).txConverg = i/casos*100;
resultadoFinal(1).txIndicador = sum(I)/N_EC2*100;