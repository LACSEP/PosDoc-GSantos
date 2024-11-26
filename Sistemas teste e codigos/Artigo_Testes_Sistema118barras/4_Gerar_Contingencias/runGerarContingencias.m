clc; clear;

%% Dados do sistema teste
addpath('.\Sistemas')

system = 'case118bus_modificado.cdf';
[nbus,NA,~,pq,pv,ref,gen,Ybus,Vm,Va,Pg0,Qg,Pd0,Qd0,Qmax,Qmin,Vmax,Vmin,Sb,names,buses,dbranch] = systemData(system);

%Barras pv, pq e geradores (inclui slack)
gen = sort(gen); nger = length(gen);

%Ramos do sistema
nBranch = length(dbranch);

%Eolicas
wind = find(Pd0<0);  nWind = length(wind); %Buses with wind farm

exclusaoLinhas = [];
for i = 1 : 1 : nger
    ilhamento(i) = 2;
    for j = 1 : 1 : 186 %186 = nBranch - nWind
        if dbranch(j).iBus == gen(i) || dbranch(j).jBus == gen(i)
            if ilhamento(i) > 0
                ilhamento(i) = ilhamento(i)-1;
            end
        end
    end
    if ilhamento(i) == 1
        j = 0;
        while true
            j = j + 1;
            if dbranch(j).iBus == gen(i) || dbranch(j).jBus == gen(i)
                exclusaoLinhas = [exclusaoLinhas j];
                break;
            end
        end
    end
end

contigencias = [1 : 1 : (nBranch-nWind)];
contigencias(exclusaoLinhas) = [];
nContigencias = length(contigencias);

fullpathname = strcat('.\Sistemas\',system);

fileID = fopen(fullpathname,'rt+');

i = 0;
while true
    i = i+1;
    tline = fgetl(fileID);
    A{i} = tline;
    if ischar(tline) == 0
        break;
    end
end
fclose(fileID);
lDataFile = length(A);

for i = 1 : 1 : lDataFile
    if contains(string(A{i}),'BRANCH DATA FOLLOWS') == 1
        caso0  = i;
    end
end

A{1} = A{1}(1:38);
for i = 1 : 1 : nContigencias
    k = 0;
    caso = caso0 + contigencias(i);
    fullpathname2 = strcat('.\Sistemas\Contingencias\Linha',num2str(contigencias(i)),'_',num2str(dbranch(contigencias(i)).iBus),'_',num2str(dbranch(contigencias(i)).jBus),'.cdf');
    B = {};
    for j = 1 : 1 : lDataFile        
        if j ~= caso
           % if j == 1 
                k = k + 1;
            %    B{k} = strcat(string(A{k}),'- Contigencia na linha ',num2str(dbranch(contigencias(i)).iBus),'-',num2str(dbranch(contigencias(i)).jBus));
            %else
             %   k = k + 1;
                B{k} = A{j};
            %end
        end
    end
    %B{1} = strcat(string(A{1}),'- Contigencia na linha ',num2str(dbranch(contigencias(i)).iBus),'-',num2str(dbranch(contigencias(i)).jBus));
    fout = fopen(fullpathname2, 'wt+');
    for j = 1:numel(B)
        if j == 1
            aux = strcat(string(A{1}),'- Contigencia na linha',{' '},num2str(dbranch(contigencias(i)).iBus),'-',num2str(dbranch(contigencias(i)).jBus));
            fprintf(fout,'%s\n', aux);
        elseif B{j} == -1
            break
        else
            fprintf(fout,'%s\n', B{j});
        end
    end
    fclose(fout);
end



