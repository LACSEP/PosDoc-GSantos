clc; clear;

warning off;

addpath('.\Sistemas\Contingencias\')

addpath('.\Sistemas')
system = 'case118bus_modificado.cdf';
[~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,dbranch] = systemData(system);

for linhaContingencia = 108
    %% Sistema modificado com o desligamento da linha
    system = strcat('Linha',num2str(linhaContingencia),'_',num2str(dbranch(linhaContingencia).iBus),'_',num2str(dbranch(linhaContingencia).jBus),'.cdf');

    %% ---- Leitura do sistema teste e obtenção da matriz de correlação ---- %%
    run_initialization;
    Vm_old = Vm;
    Va_old = Va;
    Pd0total = sum(Pd0(loadBuses));

    %% ---- Método da entropia cruzada ---- %%

    fprintf('----------------------MÉTODO DA EC----------------------\n')
    %% Teste preliminar - obtenção da g(x)

    %Parâmetros da entropia
    N_EC = 1e3; %no. de amostras
    N_el = ceil(0.1*N_EC); %no. de elite
    itMax_EC = 100;

    %Alvo
    percViol = 0.04;

    tic; rng(1);

    it_EC = 0; casos = 0;

    while true

        it_EC = it_EC + 1;

        % --- Amostragem - Na primeira iteração, o sorteio é feito em X e nas demais em Z
        for i = 1 : 1 : N_EC

            while true %Os sorteios sao repetidos ate o criterio de avaliacao ser satisfeito

                %Sorteio

                casos = casos + 1;

                while true

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
                Pg0(1:end ~= ref) = Pg0_old(1:end ~= ref) + Pg0_old(1:end ~= ref)/Pgt_old*(2*7+sum(Pd0(wind)));
                Pg0(ref) = Pg0_old(ref);
                Qg = zeros(nbus,1);

                %2) Cálculo da margem de carrregamento se o FP convergir
                %Fluxo de carga
                Vm = Vm_old; Va = Va_old;
                [Vm, Va, Pg, Qg, itPF] = runPowerFlow(Vm,Va,Vr,Pd0,Qd0,Pg0,Qg0,nbus,npv,ref,gen,nger,Qg0(gen),Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,Ybus,tepa,tepr,itPFmax);
                if itPF < itPFmax
                    %Método Direto
                    [tload, dPtotal, itPoC] = runDirectMethod(Vm,Va,Vr,Pd0,Pg,Qg,nbus,npv,ref,gen,nger,Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,dA_dQg,dC_dVa,Ybus,bpl,tan0,fpart,tol,itPoCmax);
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

        if it_EC == 1
            fprintf('Resultado do método de EC: \n');
        end

        fprintf('Iteracao %d da EC com %.2f%% dos casos convergidos: tr = %.4f \n', it_EC, i/casos*100, tr)

        if  tr < percViol || it_EC == itMax_EC
            fEC(:,1) = S_EC(ix,1);
            fEC(:,2) = S_EC(ix,2);
            fEC(:,3) = S_EC(ix,3);
            break;
        end

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

    if it_EC > 1

        N_EC2 = 1e4;

        I = []; W = []; p = [];

        rng(1); d = 10; casos = 0;

        for i = 1 : 1 : N_EC2

            while true

                casos = casos + 1;

                %Sorteio
                while true
                    X2(i,:) = mvnrnd(muX,sigmaX,1);
                    if min(X2(i,:)) >= 0 %Garante que carga continua consumindo e eólica gerando
                        break;
                    end
                end

                % Avaliação

                %1) Atualização dos valores de carga e geração
                bpl(loadBuses) = X2(i,1:nLoad);
                Pd0(wind) = -1*X2(i,nW1:end);
                Pd0total = sum(bpl);
                Pg0(1:end ~= ref) = Pg0_old(1:end ~= ref) + Pg0_old(1:end ~= ref)/Pgt_old*(2*7+sum(Pd0(wind)));
                Pg0(ref) = Pg0_old(ref);
                Qg = zeros(nbus,1);

                %2) Cálculo da margem de carrregamento se o FP convergir
                %Fluxo de carga
                Vm = Vm_old; Va = Va_old;
                [Vm, Va, Pg, Qg, itPF] = runPowerFlow(Vm,Va,Vr,Pd0,Qd0,Pg0,Qg0,nbus,npv,ref,gen,nger,Qg0(gen),Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,Ybus,tepa,tepr,itPFmax);
                if itPF < itPFmax
                    %Método Direto
                    [tload, dPtotal, itPoC] = runDirectMethod(Vm,Va,Vr,Pd0,Pg,Qg,nbus,npv,ref,gen,nger,Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,dA_dQg,dC_dVa,Ybus,bpl,tan0,fpart,tol,itPoCmax);
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

        erroCE = 2.5758*std(p)/sqrt(N_EC2)/mean(p);   %2.5758 = icdf('norm',(1-0.99)/2,0,1)

    else

        N_MCS = 1e4;

        I = []; W = []; p = [];

        rng(1); d = 10; casos = 0;

        for i = 1 : 1 : N_MCS

            while true

                casos = casos + 1;

                %Sorteio
                while true
                    X2(i,:) = mvnrnd(mu,Sigma,1);
                    if min(X2(i,100:106)) >= 0 %Garante que carga continua consumindo e eólica gerando
                        break;
                    end
                end

                % Avaliação

                %1) Atualização dos valores de carga e geração
                bpl(loadBuses) = X2(i,1:nLoad);
                Pd0(wind) = -1*X2(i,nW1:end);
                Pd0total = sum(bpl);
                Pg0(1:end ~= ref) = Pg0_old(1:end ~= ref) + Pg0_old(1:end ~= ref)/Pgt_old*(2*7+sum(Pd0(wind)));
                Pg0(ref) = Pg0_old(ref);
                Qg = zeros(nbus,1);

                %2) Cálculo da margem de carrregamento se o FP convergir
                %Fluxo de carga
                Vm = Vm_old; Va = Va_old;
                [Vm, Va, Pg, Qg, itPF] = runPowerFlow(Vm,Va,Vr,Pd0,Qd0,Pg0,Qg0,nbus,npv,ref,gen,nger,Qg0(gen),Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,Ybus,tepa,tepr,itPFmax);
                if itPF < itPFmax
                    %Método Direto
                    [tload, dPtotal, itPoC] = runDirectMethod(Vm,Va,Vr,Pd0,Pg,Qg,nbus,npv,ref,gen,nger,Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,dA_dQg,dC_dVa,Ybus,bpl,tan0,fpart,tol,itPoCmax);
                    %Verifica convergência do método direto
                    if  itPoC <= itPoCmax
                        Pd = Pd0 + tload*bpl;
                        %Verifica se as cargas continuam maiores do que 0
                        if min(Pd(loadBuses)) >= 0
                            Pd = Pd0 + tload*bpl;
                            S_MC(i,1) = dPtotal/Pd0total;
                            S_MC(i,2) = tload;
                            S_MC(i,3) = Pd0total;
                            break;
                        end
                    end
                end
            end

            if S_MC(i,1) < percViol
                I(i) = 1;                
            else
                I(i) = 0;                
            end

            if sum(I) == 0
                p(i) = 0;
            else
                p(i) = 1/i*sum(I);
            end


            if rem(i,d) == 0
                if d == 10
                    fprintf('----------------------Resultados parciais do método de Monte Carlo----------------------\n')
                end
                fprintf('Para %d iterações: %.4f%% com %.2f dos casos convergidos \n', d, p(end)*100, i/casos*100);
                d = d * 10;
            end
        end
        erroMCS = 2.5758*std(p)/sqrt(N_MCS)/mean(p);   %2.5758 = icdf('norm',(1-0.99)/2,0,1)
    end
end
%%
toc;
iteste = 1;
resultadoFinal(iteste).linha = linhaContingencia;
resultadoFinal(iteste).BarraDe = dbranch(linhaContingencia).iBus;
resultadoFinal(iteste).BarraPara = dbranch(linhaContingencia).jBus;
resultadoFinal(iteste).SMC_preAmost = N_MCS;
resultadoFinal(iteste).SMC_4porc = p(end)*100;
resultadoFinal(iteste).tempo = toc;
resultadoFinal(iteste).erro = erroMCS;
resultadoFinal(iteste).txConvergencia = i/casos*100;

