%% ---- Método da entropia cruzada ---- %%

fprintf('----------------------MÉTODO DA EC----------------------\n')
%% Teste preliminar - obtenção da g(x)
tic;

it_EC = 0;

Pd0total = sum(Pd0(loadBuses));

X = []; U = []; Z = []; S_EC = zeros(N_EC,3);

rng(2);

while true

    it_EC = it_EC + 1;

    casos = 0;
    % --- Amostragem - Na primeira iteração, o sorteio é feito em X e nas demais em Z
    for i = 1 : 1 : N_EC

        while true %Os sorteios sao repetidos ate o criterio de avaliacao ser satisfeito

            %Sorteio

            %while true

            casos = casos + 1;

            if it_EC == 1 %Sorteio em X
                X(i,:) = mvnrnd(mu,Sigma,1);

            else %Sorteio em Z
                Z(i,:) = mvnrnd(muZ,sigmaZ,1);
                for j = 1 : 1 : nLoad
                    U(i,j) = normcdf(Z(i,j));
                end
                for j = 1 : 1 : nLoad
                    Z(i,j) = norminv(U(i,j),0,1);
                end
                for j = 1 : 1 : nLoad
                    X(i,j) = Z(i,j) * sigX(j) + muX(j);
                end

            end

            %             if min(X(i,:)) >= 0 %Garante que carga continua consumindo e eólica gerando
            %                 break;
            %             end

            %         end

            % Avaliação

            %1) Atualização dos valores de carga e geração
            bpl(loadBuses) = X(i,:);

            %2) Cálculo da margem de carrregamento se o FP convergir
            %Fluxo de carga
            [Vm, Va, Pg, Qg, itPF] = runPowerFlow(Vm0,Va0,Vr,Pd0,Qd0,Pg0,Qg0,nbus,npv,ref,gen,nger,Qg0(gen),Qmax,Qmin,dP_dQg,dQ_dQg,dy_dVa,Ybus,tepa,tepr,itPFmax);
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
