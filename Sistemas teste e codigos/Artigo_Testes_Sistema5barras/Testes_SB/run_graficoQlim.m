clc; clear;

warning off;

% casebase -> CARTAO NO FORMATO MATPOWER!
casobase = case5;
[nbus,narea,area,pq,pv,ref,gen,Ybus,Vm,Va,Pg0,Qg,Pd0,Qd0,Qmax2,Qmin2,Vmax,Vmin,Sb,names,buses] = data(casobase);
Qmax = zeros(nbus,1); Qmax(gen) = Qmax2;
Qmin = zeros(nbus,1); Qmin(gen) = Qmin2;

%Ordering and size of data vectors
npv = length(pv); npq = length(pq);
gen = sort(gen); nger = length(gen);

%Reactive power limit of SGs
Qmax = zeros(nbus,1); Qmax(gen) = Qmax2;
Qmin = zeros(nbus,1); Qmin(gen) = Qmin2;
if Qmax(ref) == 0   Qmax(ref) = 999.99;     end
if Qmin(ref) == 0   Qmin(ref) = -999.99;    end

%Initial Values
Vr = Vm;
Qspc = Qg(gen);

%Control Constants
sigk = 1e8; sigv = 1e-6; sigq = 1e-6;
tepa = 1e-8; tepr = 1e-8; qlst = 0.04;

%Matrices without changes
dP_dQg = zeros(nbus-1,nger);
dQ_dQg = zeros(nbus,nger);
for i = 1 : 1 : nger
    j = gen(i);
    dQ_dQg(j,i) = -1;
end
dy_dVa = zeros(nger,nbus-1);
dA_dQg = zeros(nbus-1,nger);
dC_dVa = zeros(nger,nbus-1);

for casoQmax = 6
    Qmax(3) = casoQmax;

    for casofp = 0 : 1 : 10

        %% ---- Base case power flow ---- %%
        powerFlow;
        
        %Pg and Qspc update
        Pg0(gen) = (real(Scalc(gen))-Pd0(gen));
        Qspc = Qg(gen);

        %% ---- Direct method parameters ---- %%
        %Identification of loaded and unloaded buses
        Pd0_total = sum(Pd0);
        Sd0 = sqrt(Pd0.^2+Qd0.^2);
        pos = find(Sd0==0); npos = length(pos); %Unloaded buses
        pos2 = find(Sd0~=0); npos2 = length(pos2); %Loaded buses
        %Dispatch of SGs
        fpart = zeros(nbus,1);
        % fpart(gen) = Pg0(gen)/sum(Pg0);
        fpart(1) = casofp/10;
        fpart(3) = 1-fpart(1);

        %Load growth factor
        bp0 = zeros(nbus,1);
        bp0(pos2) = Pd0(pos2)/sum(Pd0);
        alphaq = Qd0./Pd0; alphaq(pos) = 0; %Qd0 é aumentado de forma a manter fp cte;
        %Initial Values
        t = 0; %Load Growth Parameter
        wP = ones(nbus,1); wP(ref) = 0;
        wQ = ones(nbus,1);
        wG = ones(nger,1);
        tol = 1e-7;  itMax = 200; %Tolerance (tolerance and maximum number of iterations)
        casof = 0; casog = 0; bp = zeros(nbus,1);
        for i = 0 : 0.1 : 1
            for j = 0 : 0.1 : 1
                for k = 0 : 0.1 : 1
                    if i ~= 0 || j ~= 0 || k ~= 0
                        bp(pos2(1)) = bp0(pos2(1))*i;
                        bp(pos2(2)) = bp0(pos2(2))*j;
                        bp(pos2(3)) = bp0(pos2(3))*k;
                        [tload,dPtotal,it] = run_DirectMethodQlim(t, wP, wQ, wG, tol, gen, nger, pq, npq, ref, pv, npv,Ybus,Vm,Va,Pg0,Qg,Qmin,Qmax,Pd0,nbus,itMax,bp,alphaq,fpart);
                        if  it <= itMax
                            Pd = Pd0 + tload*bp;
                            if min(Pd) >= 0
                                casof = casof + 1;
                                f(casof,1) = dPtotal/Pd0_total;
                                for p = 1 : 1 : npos2
                                    f(casof,p+1) = Pd0(pos2(p)) + bp(pos2(p))*tload;
                                end
                            end
                        end
                        [tload,dPtotal,it] = run_DirectMethod(t, wP, wQ, wG, tol, gen, nger, pq, npq, ref, pv, npv,Ybus,Vm,Va,Pg0,Qg,Qmin,Qmax,Pd0,nbus,itMax,bp,alphaq,fpart);
                        if  it <= itMax
                            Pd = Pd0 + tload*bp;
                            if min(Pd) >= 0
                                casog = casog + 1;
                                g(casog,1) = dPtotal/Pd0_total;
                                for p = 1 : 1 : npos2
                                    g(casog,p+1) = Pd0(pos2(p)) + bp(pos2(p))*tload;
                                end
                            end
                        end
                    end
                end
            end
        end

        % Create figure
        figure('PaperSize',[30 20],'WindowState','maximized');

        % Create axes
        axes1 = axes;
        hold(axes1,'on');
        grid(axes1,'on');
        axis(axes1,'tight');
        hold(axes1,'off');


        xlin = linspace(min(g(:,2)),max(g(:,2)),length(g));
        ylin = linspace(min(g(:,3)),max(g(:,3)),length(g));
        [X,Y] = meshgrid(xlin,ylin);
        Z = griddata(g(:,2),g(:,3),g(:,4),X,Y,'natural');
        mesh(X,Y,Z,'EdgeColor',[0.149019607843137 0.149019607843137 0.149019607843137])
        hold on

        xlin = linspace(min(f(:,2)),max(f(:,2)),length(f));
        ylin = linspace(min(f(:,3)),max(f(:,3)),length(f));
        [X,Y] = meshgrid(xlin,ylin);
        Z = griddata(f(:,2),f(:,3),f(:,4),X,Y,'natural');
        mesh(X,Y,Z,'EdgeColor',[0.501960784313725 0.501960784313725 0.501960784313725]);

        axis tight;

        % Create zlabel
        zlabel('PL_5 (p.u.)','FontName','Times');

        % Create ylabel
        ylabel('PL_4 (p.u.)','FontName','Times');

        % Create xlabel
        xlabel('PL_2 (p.u.)','FontName','Times');

        legend('Without Q-limits', 'With Q-limits')

        %view(axes1,[20.6765141803106 10.6897164040909]); %Qlim = 2
        %view(axes1,[37.3015143002093 7.86094615233248]); %Qlim = 3
        %view(axes1,[31.4390139513203 5.79491294953289]); %Qlim = 4
        %view(axes1,[-5.48598709459683 27.888433426775]); %Qlim = 5
        view(axes1,[-38.2109860377043 13.3635000205435]); %Qlim = 6


        title(['SB para Qmax = ',num2str(Qmax(3)),' p.u.'],['FP = [GS_1: ',num2str(fpart(1)),'- GS_3:',num2str(fpart(3)),']']);

        % Set the remaining axes properties
        set(axes1,'FontName','Times','FontSize',30);

        pause(3);

        % Create legend
        legend1 = legend(axes1,'show');
        set(legend1,...
            'Position',[0.712217017358608 0.729712044443762 0.155208329499389 0.0905759136714235]);
        print(['SB_Qmax_',num2str(Qmax(3)),'_CasoFPart_',num2str(casofp),'.pdf'],'-dpdf','-r200')

        %savefig(['SB_Qmax_',num2str(Qmax(3)),'_CasoFPart_',num2str(casofp),'.fig']);

        pause(3);

        clf; close all;

    end

end

