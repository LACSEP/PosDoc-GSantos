%% Step 2 - MCS from CE results

it_EC2 = 0;

I = zeros(1,N_EC2); W = zeros(1,N_EC2); p = zeros(1,N_EC2); 

rng(2); 

d = 10; casos = 0;

for i = 1 : 1 : N_EC2

    while true

        %Sorteio
       % while true
            casos = casos + 1;
            X2(i,:) = mvnrnd(muX,sigmaX,1);
         %   if min(X2(i,:)) >= 0 %Garante que carga continua consumindo e eólica gerando
            %    break;
         %   end
     %   end

        % Avaliação

        %1) Atualização dos valores de carga e geração
        bpl(loadBuses) = X2(i,1:nLoad);

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
        f(i,1) = mvnpdf(X2(i,:),mu,Sigma);
        g(i,1) = mvnpdf(X2(i,:),muX,sigmaX);
        W(i) = f(i,1)/g(i,1);
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
        fprintf('Para %d iterações: %.4f%% com %.2f dos casos convergidos \n', d, p(i)*100, i/casos*100);
        d = d * 10;
    end
end
