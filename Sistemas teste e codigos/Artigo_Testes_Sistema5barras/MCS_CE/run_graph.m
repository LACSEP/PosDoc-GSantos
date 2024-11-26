clc; clear;


for i = 3


    %i = 1 - Gráfico sem correlação

    if i == 1

        X1 = load('.\Resultados\SMC_rho0.0_3e3.mat');
        X1 = X1.S_MC(:,1);

        X2 = load('.\Resultados\SMC_rho0.0_1e4.mat');
        X2 = X2.S_MC(:,1);

        X3 = load('.\Resultados\SMC_rho0.0.mat');
        X3 = X3.S_MC(:,1);
        
        X4 = load('.\Resultados\SEC_rho0.0.mat');
        X4 = X4.S2_EC(:,1);
        W = load('.\Resultados\WEC_rho0.0.mat');
        W = W.W';

        minX = min(min(min(X1),min(X2)),min(min(X3),min(X4)));
        maxX = 0.07;

        nbin = 5;
        %edges = linspace(minX,0.07,nbin); 
        edges = [minX, 0.049, 0.056, 0.063, 0.07];
        
        
        for j = 2 : 1 : nbin
            bin_mask1 = X1 >= edges(j-1) & X1 < edges(j);
            bin_mask2 = X2 >= edges(j-1) & X2 < edges(j);
            bin_mask3 = X3 >= edges(j-1) & X3 < edges(j);
            bin_mask4 = X4 >= edges(j-1) & X4 < edges(j);
            bin_counts1(j-1) = sum(bin_mask1)/3e3*100; 
            bin_counts2(j-1) = sum(bin_mask2)/1e4*100; 
            bin_counts3(j-1) = sum(bin_mask3)/1e6*100;
            bin_counts4(j-1) = sum(W(bin_mask4))/3e3*100;
            
            if j == 2
                dleg{j-1} = strcat('<',num2str(edges(j)*100,'%.1f'),'%');
            else
                dleg{j-1} = strcat(num2str(edges(j-1)*100,'%.1f'),'-',num2str(edges(j)*100,'%.1f'),'%');
            end
        end            

        y = [bin_counts1(1),bin_counts1(2),bin_counts1(3),bin_counts1(4);bin_counts2(1),bin_counts2(2),bin_counts2(3),bin_counts2(4);bin_counts3(1),bin_counts3(2),bin_counts3(3),bin_counts3(4);bin_counts4(1),bin_counts4(2),bin_counts4(3),bin_counts4(4)];

        b = bar(y,'stacked');        

        clrs = [0.2 0.2 0.2; 0.45 0.45 0.45; 0.65 0.65 0.65; 0.9 0.9 0.9];
        set(b,{'FaceColor'},{clrs(1,:),clrs(2,:),clrs(3,:),clrs(4,:)}.')

        grid;

        legend(dleg)

        xlabel('(a)');

        ylabel('Probability of violation (%)');

        %Set new xtick and xticklabels, rotate by 90deg.
        ax = b.Parent; % axis handle, if you don't have it already
        ax.XTick = [1 2 3 4];      
        set(gca(),'XTickLabel',{sprintf('MCS w/ 3\\cdot10^3 \\newlinesamples') sprintf('MCS w/ 10^4 \\newlinesamples') sprintf('MCS w/ 10^6 \\newlinesamples') sprintf('CE w/ 0.75/3\\cdot10^3\\newlinesamples')});        
        ax.XTickLabelRotation = 45;      

        set(gca,'FontSize',33,'FontName','Times')
        set(gcf,'Paperunits','inches','PaperPosition',[0 0 15 10],'PaperSize',[15 10])
        pause(1);
        print('Grafico1_2.pdf','-dpdf','-r400')

    elseif i == 2

        X1 = load('.\Resultados\SMC_rho0.4_3e3.mat');
        X1 = X1.S_MC(:,1);

        X2 = load('.\Resultados\SMC_rho0.4_1e4.mat');
        X2 = X2.S_MC(:,1);

        X3 = load('.\Resultados\SMC_rho0.4.mat');
        X3 = X3.S_MC(:,1);
        
        X4 = load('.\Resultados\SEC_rho0.4.mat');
        X4 = X4.S2_EC(:,1);
        W = load('.\Resultados\WEC_rho0.4.mat');
        W = W.W';

        minX = min(min(min(X1),min(X2)),min(min(X3),min(X4)));
        maxX = 0.07;

        nbin = 5;
        %edges = linspace(minX,0.07,nbin); 
        edges = [minX, 0.049, 0.056, 0.063, 0.07];
        
        
        for j = 2 : 1 : nbin
            bin_mask1 = X1 >= edges(j-1) & X1 < edges(j);
            bin_mask2 = X2 >= edges(j-1) & X2 < edges(j);
            bin_mask3 = X3 >= edges(j-1) & X3 < edges(j);
            bin_mask4 = X4 >= edges(j-1) & X4 < edges(j);
            bin_counts1(j-1) = sum(bin_mask1)/3e3*100; 
            bin_counts2(j-1) = sum(bin_mask2)/1e4*100; 
            bin_counts3(j-1) = sum(bin_mask3)/1e6*100;
            bin_counts4(j-1) = sum(W(bin_mask4))/3e3*100;
            
            if j == 2
                dleg{j-1} = strcat('<',num2str(edges(j)*100,'%.1f'),'%');
            else
                dleg{j-1} = strcat(num2str(edges(j-1)*100,'%.1f'),'-',num2str(edges(j)*100,'%.1f'),'%');
            end
        end            

        y = [bin_counts1(1),bin_counts1(2),bin_counts1(3),bin_counts1(4);bin_counts2(1),bin_counts2(2),bin_counts2(3),bin_counts2(4);bin_counts3(1),bin_counts3(2),bin_counts3(3),bin_counts3(4);bin_counts4(1),bin_counts4(2),bin_counts4(3),bin_counts4(4)];

        b = bar(y,'stacked');        

        clrs = [0.2 0.2 0.2; 0.45 0.45 0.45; 0.65 0.65 0.65; 0.9 0.9 0.9];
        set(b,{'FaceColor'},{clrs(1,:),clrs(2,:),clrs(3,:),clrs(4,:)}.')

        grid;

       % legend(dleg)

        xlabel('(b)');

        ylabel('Probability of violation (%)');        
        
        %Set new xtick and xticklabels, rotate by 90deg.
        ax = b.Parent; % axis handle, if you don't have it already
        ax.XTick = [1 2 3 4];      
        set(gca(),'XTickLabel',{sprintf('MCS w/ 3\\cdot10^3 \\newlinesamples') sprintf('MCS w/ 10^4 \\newlinesamples') sprintf('MCS w/ 10^6 \\newlinesamples') sprintf('CE w/ 0.75/3\\cdot10^3\\newlinesamples')});        
        ax.XTickLabelRotation = 45;      

        set(gca,'FontSize',33,'FontName','Times')
        set(gcf,'Paperunits','inches','PaperPosition',[0 0 15 10],'PaperSize',[15 10])
        pause(1);
        print('Grafico2_2.pdf','-dpdf','-r400')

    else

        X1 = load('.\Resultados\SMC_rho0.8_3e3.mat');
        X1 = X1.S_MC(:,1);

        X2 = load('.\Resultados\SMC_rho0.8_1e4.mat');
        X2 = X2.S_MC(:,1);

        X3 = load('.\Resultados\SMC_rho0.8.mat');
        X3 = X3.S_MC(:,1);
        
        X4 = load('.\Resultados\SEC_rho0.8.mat');
        X4 = X4.S2_EC(:,1);
        W = load('.\Resultados\WEC_rho0.8.mat');
        W = W.W';

        minX = min(min(min(X1),min(X2)),min(min(X3),min(X4)));
        maxX = 0.07;

        nbin = 5;
        %edges = linspace(minX,0.07,nbin); 
        edges = [minX, 0.049, 0.056, 0.063, 0.07];
        
        
        for j = 2 : 1 : nbin
            bin_mask1 = X1 >= edges(j-1) & X1 < edges(j);
            bin_mask2 = X2 >= edges(j-1) & X2 < edges(j);
            bin_mask3 = X3 >= edges(j-1) & X3 < edges(j);
            bin_mask4 = X4 >= edges(j-1) & X4 < edges(j);
            bin_counts1(j-1) = sum(bin_mask1)/3e3*100; 
            bin_counts2(j-1) = sum(bin_mask2)/1e4*100; 
            bin_counts3(j-1) = sum(bin_mask3)/1e6*100;
            bin_counts4(j-1) = sum(W(bin_mask4))/3e3*100;
            
            if j == 2
                dleg{j-1} = strcat('<',num2str(edges(j)*100,'%.1f'),'%');
            else
                dleg{j-1} = strcat(num2str(edges(j-1)*100,'%.1f'),'-',num2str(edges(j)*100,'%.1f'),'%');
            end
        end            

        y = [bin_counts1(1),bin_counts1(2),bin_counts1(3),bin_counts1(4);bin_counts2(1),bin_counts2(2),bin_counts2(3),bin_counts2(4);bin_counts3(1),bin_counts3(2),bin_counts3(3),bin_counts3(4);bin_counts4(1),bin_counts4(2),bin_counts4(3),bin_counts4(4)];

        b = bar(y,'stacked');        

        clrs = [0.2 0.2 0.2; 0.45 0.45 0.45; 0.65 0.65 0.65; 0.9 0.9 0.9];
        set(b,{'FaceColor'},{clrs(1,:),clrs(2,:),clrs(3,:),clrs(4,:)}.')

        grid;

       % legend(dleg)

        xlabel('(c)');

        ylabel('Probability of violation (%)');

        %Set new xtick and xticklabels, rotate by 90deg.
        ax = b.Parent; % axis handle, if you don't have it already
        ax.XTick = [1 2 3 4];      
        set(gca(),'XTickLabel',{sprintf('MCS w/ 3\\cdot10^3 \\newlinesamples') sprintf('MCS w/ 10^4 \\newlinesamples') sprintf('MCS w/ 10^6 \\newlinesamples') sprintf('CE w/ 0.75/3\\cdot10^3\\newlinesamples')});        
        ax.XTickLabelRotation = 45;      
        
        set(gca,'FontSize',33,'FontName','Times')
        set(gcf,'Paperunits','inches','PaperPosition',[0 0 15 10],'PaperSize',[15 10])
        pause(1);
        print('Grafico3_2.pdf','-dpdf','-r400')

    end


end



