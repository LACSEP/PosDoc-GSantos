clc; clear;

warning off;

addpath('.\Sistemas')

%Sistema
system = 'case118bus_modificado.cdf';

%% ---- Leitura do sistema teste e obtenção da matriz de correlação ---- %%
run_initialization;

% %% ---- Mapa de calor da matriz de correlação ---- %%
% heatmap(corrSigma);
% xlabel('Loads 1-99       WFs 1-7');
% ylabel('Loads 1-99       WFs 1-7');% 
% set(gca,'FontSize',20,'FontName','Times');
% set(gcf,'Paperunits','inches','PaperPosition',[0 0 20 20],'PaperSize',[20 20])
% xtick = {};
% pause(1);
% print('Heatmap_corrSigma.pdf','-dpdf','-r400');

%% ---- Mapa de calor da matriz de correlação entre eolicas ---- %%
heatmap(corrSigma(100:106, 100:106), 'XData', ["8","32","42","55","76","92","105"], 'YData',["8","32","42","55","76","92","105"]);
xlabel('Wind-generation buses');
ylabel('Wind-generation buses');% 

set(gca,'FontSize',60,'FontName','Times');
set(gcf,'Paperunits','inches','PaperPosition',[0 0 20 20],'PaperSize',[20 20])
pause(1);
colormap(flip(gray))

%print('Heatmap_corrSigmaWF.pdf','-dpdf','-r400');

