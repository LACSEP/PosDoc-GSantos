clc; clear;

%% Dados do sistema teste
addpath('.\Sistemas')
system = 'case118bus_modificado.cdf';
[nbus,~,~,~,~,~,gen,~,~,~,~,~,Pd0,~,~,~,~,~,~,~,~,dbranch] = systemData(system);

%Barras pv, pq e geradores (inclui slack)
gen = sort(gen); nger = length(gen);

%Ramos do sistema
nBranch = length(dbranch);

%Eolicas
wind = find(Pd0<0);  nWind = length(wind); %Buses with wind farm

%% Lista de contingencias
exclusaoLinhas = [];
n = nBranch - nWind;
for i = 1 : 1 : nger
    ilhamento(i) = 2;
    for j = 1 : 1 : n %186 = nBranch - nWind
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

exclusaoLinhas = [exclusaoLinhas];

contingencias = [1 : 1 : (nBranch-nWind)];
contingencias(exclusaoLinhas) = [];
nContigencias = length(contingencias);

[filepath,name,ext] = fileparts(system);

fullpathname = strcat('.\Sistemas\',name,'.pwf'); %Arquivo .pwf gerado para as contingencias

fileID = fopen(fullpathname,'rt+');

i = 0;
while true
    i = i+1;
    tline = fgetl(fileID);
    A{i} = tline;   
    if contains(string(A{i}),'TITU') == 1
        iTitu = i + 1;
    end
    if ischar(tline) == 0
        break;
    end
end
fclose(fileID);


for i = 1 : 1 : nContigencias
    B = A; d = 0;
    B{iTitu} = strcat('Sistema de',32,num2str(nbus-nWind),' barras - Contigencia na linha', 32,num2str(contingencias(i)),32,'- da barra',32,num2str(dbranch(contingencias(i)).iBus),32,'para barra',32,num2str(dbranch(contingencias(i)).jBus),32);
    j = 151 + contingencias(i);
    B{j}(18) = 'D';  

    nomeArquivo = strcat('G:\Meu Drive\SELECAO\Artigo_Testes_Sistema118barras\5_ExecutarPoC_ANAREDE\Sistemas\ContingenciasANAREDE\Linha',num2str(contingencias(i)),'_',num2str(dbranch(contingencias(i)).iBus),'_',num2str(dbranch(contingencias(i)).jBus),'.pwf');
    
    fout = fopen(nomeArquivo, 'wt+');

    for k = 1:numel(B)-1
        if B{k} == -1
            fprintf(fout,'%s', B{k});
            break
        else
            fprintf(fout,'%s\n', B{k});
        end
    end

    fclose(fout);


    disp('Executando ANAREDE...')


    filepath = ['"' nomeArquivo '"']

    eval(['!start /MIN cmd /c C:\\CEPEL\\Anarede\\V110702\\ANAREDE.exe',char(32),filepath]);

    pause(240)
    
    eval(['!start /min cmd /q /c taskkill /f /im ANAREDE.exe'])
    
    delete('*.sav')

    fullpathname = 'G:\Meu Drive\SELECAO\Artigo_Testes_Sistema118barras\5_ExecutarPoC_ANAREDE\Sistemas\ContingenciasANAREDE\Relat.out';

    fileID = fopen(fullpathname,'rt+');
    k = 0;
    while true
        k = k+1;
        tline = fgetl(fileID);
        C{k} = tline;        
        if ischar(tline) == 0
            break;
        end
        if contains(C{k},'Carregamento do sistema =')
            d = k;
            resultado(i).No = i;
            resultado(i).Linha = contingencias(i);
            resultado(i).DeBarra = dbranch(contingencias(i)).iBus;
            resultado(i).ParaBarra = dbranch(contingencias(i)).jBus;
            resultado(i).Incremento = str2double(regexp(C{k},'[\d.]+','match'));
        end
    end
    fclose(fileID);
    delete('*.out');
    pause(10)
    d = d + 6;
    for j = 1 : 1 : 10
        resultado(i).Relat2{j} = C{d+j};
    end   
end

%%
j = 0; 
for i = 1 : 1 : nContigencias
    j = j + 1;
    D{j} = strcat('----------CASO',32,num2str(i),'----------');
    j = j + 1;
    D{j} = strcat('Contingência na linha',32,num2str(contingencias(i)));
    j = j + 1;
    D{j} = strcat('Barra De:',32,num2str(dbranch(contingencias(i)).iBus));
    j = j + 1;
    D{j} = strcat('Barra Para:',32,num2str(dbranch(contingencias(i)).jBus));
    j = j + 1;
    D{j} = strcat('Incremento de carregamento do sistema:',32,num2str(resultado(i).Incremento,'%4.4f'),'%');
    j = j + 1;
    D{j} = 'X------------------X-------------------X'
    j = j + 1;
    D{j} = '      Barra                             ';
    j = j + 1;
    D{j} = ' Num.      Nome      Tensao   Variacao  ';
    j = j + 1;
    D{j} = 'X-----X------------X---------X---------X';
    for k = 1 : 1 : 10
        j = j + 1;
        D{j} = resultado(i).Relat2{k}
    end
    j = j + 1;
end

nomeArquivo = 'C:\Users\santo\Meu Drive\POS-DOUTORADO\ArquivosFinaisFinais\Testes_Sistema118barras\Sistemas\ContingenciasANAREDE\Resultados.txt';
    
fout = fopen(nomeArquivo, 'wt+');

for k = 1:numel(D)-1
    if D{k} == -1
        fprintf(fout,'%s', D{k});
        break
    else
        fprintf(fout,'%s\n', D{k});
    end
end

fclose(fout);


