clc; clear;

system = '107barras';
%system = 'NETS_NYPS_vb';

fullpathname = strcat('.\Sistemas\',system,'.pwf');

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

dcte.base = {};
dcte.dase = {};

n = length(A);
i = 0;
while true
    i = i + 1;
    if strfind(A{i},'TITU') == 1 %Branch Data
        while true
            i = i + 1;
            if strcmp(strtrim(A{i}(1)),'(') == 0
                titu = strtrim(A{i}(1:end));
                break;
            end
        end
    end
    if strfind(A{i},'DCTE') == 1 %Branch Data
        while true
            i = i + 1;
            if contains(A{i}(1:5),'99999') == 1
                break;
            end
            if strcmp(strtrim(A{i}(1)),'(') == 0 %
                if isempty(strfind(A{i},'BASE')) == 0%Branch Data
                    k = strfind(A{i},'BASE');
                    dcte.base = str2num(A{i}(k+5:k+10));
                end
            end
            if strcmp(strtrim(A{i}(1)),'(') == 0 %
                if isempty(strfind(A{i},'DASE')) == 0 %Branch Data
                    k = strfind(A{i},'DASE');
                    dcte.dase = str2num(A{i}(k+5:k+10));
                end
            end
        end
    end

    if strfind(A{i},'DBAR') == 1 %Branch Data
        j = 0;
        while true
            i = i + 1;
            if contains(A{i}(1:5),'99999') == 1
                break;
            end
            if strcmp(strtrim(A{i}(1)),'(') == 0 %
                j = j + 1;
                dbar(j).numero = str2double(A{i}(1:5));
                dbar(j).operacao = A{i}(6); %(A ou 0) ou (E ou 1) ou (M ou 2)
                dbar(j).estado = A{i}(7); %L ou D
                dbar(j).tipo = str2double(A{i}(8));
                dbar(j).grupo_base_tensao = A{i}(9:10); %2 caracteres (letra ou numero)
                dbar(j).nome = A{i}(11:22);
                dbar(j).grupo_limite_tensao = str2double(A{i}(23:24)); %2 carcateres (letra ou numero)
                dbar(j).tensao = str2double(A{i}(25:28));
                dbar(j).angulo = str2double(A{i}(29:32));
                dbar(j).potencia_ativa = str2double(A{i}(33:37));
                dbar(j).potencia_reativa = str2double(A{i}(38:42));
                dbar(j).potencia_reativa_minima = str2double(A{i}(43:47));
                dbar(j).potencia_reativa_maxima = str2double(A{i}(48:52));
                dbar(j).barra_controlada = str2double(A{i}(53:58));
                dbar(j).demanda_ativa = str2double(A{i}(59:63));
                dbar(j).demanda_reativa = str2double(A{i}(64:68));
                try dbar(j).shunt_barra = str2double(A{i}(69:73)); end
                try dbar(j).area = str2double(A{i}(74:76)); end
                try dbar(j).tensao_base = str2double(A{i}(77:80))/1000; end
                try dbar(j).modo = str2double(A{i}(81)); end
                try dbar(j).agreg1 = str2double(A{i}(82:84)); end
                try dbar(j).agreg2 = str2double(A{i}(85:87)); end
                try dbar(j).agreg3 = str2double(A{i}(88:90)); end
                try dbar(j).agreg4 = str2double(A{i}(91:93)); end
                try dbar(j).agreg5 = str2double(A{i}(94:96)); end
                try dbar(j).agreg6 = str2double(A{i}(97:99)); end
                try dbar(j).agreg7 = str2double(A{i}(100:102)); end
                try dbar(j).agreg8 = str2double(A{i}(103:105)); end
                try dbar(j).agreg9 = str2double(A{i}(106:108)); end
                try dbar(j).agreg10 = str2double(A{i}(109:111)); end
            end
        end
    end

    if strfind(A{i},'DLIN') == 1 %Branch Data
        j = 0;
        while true
            i = i + 1;
            if contains(A{i}(1:5),'99999') == 1
                break;
            end
            if strcmp(strtrim(A{i}(1)),'(') == 0 %
                j = j + 1;
                dlin(j).de = str2double(A{i}(1:5));
                dlin(j).abertura_de = A{i}(6); %L ou D
                dlin(j).operacao = A{i}(8); % (A ou 0) ou (E ou 1) ou (M ou 2)
                dlin(j).abertura_para = A{i}(8); %L ou D
                dlin(j).para = str2double(A{i}(11:15));
                dlin(j).circuito = str2double(A{i}(16:17));
                dlin(j).estado = A{i}(18); %L ou D
                dlin(j).proprietario = A{i}(19); %F ou T
                dlin(j).resistencia = str2double(A{i}(21:26)); %R%
                dlin(j).reatancia = str2double( A{i}(27:32));%X%
                try dlin(j).susceptancia = str2double(A{i}(33:38)); end %MVAr
                try dlin(j).tap = str2double(A{i}(39:43)); end
                try dlin(j).tap_minimo = str2double(A{i}(44:48)); end
                try dlin(j).tap_maximo = str2double(A{i}(49:53)); end
                try dlin(j).tap_defasagem = str2double(A{i}(54:47)); end
                try dlin(j).barra_controlada = str2double(A{i}(59:64));end
                try dlin(j).capacidade_normal = str2double(A{i}(65:68)); end
                try dlin(j).capacidade_emergencial = str2double(A{i}(69:72)); end
                try dlin(j).numero_taps = str2double(A{i}(73:74)); end
                try dlin(j).capacidade_equipamento = str2double(A{i}(75:78)); end
                try dlin(j).agreg1 = str2double(A{i}(79:81)); end
                try dlin(j).agreg2 = str2double(A{i}(82:84)); end
                try dlin(j).agreg3 = str2double(A{i}(85:87)); end
                try dlin(j).agreg4 = str2double(A{i}(88:90)); end
                try dlin(j).agreg5 = str2double(A{i}(91:93)); end
                try dlin(j).agreg6 = str2double(A{i}(94:96)); end
                try dlin(j).agreg7 = str2double(A{i}(97:99)); end
                try dlin(j).agreg8 = str2double(A{i}(100:102)); end
                try dlin(j).agreg9 = str2double(A{i}(103:105)); end
                try dlin(j).agreg10 = str2double(A{i}(106:108)); end
            end
        end
    end
    if strfind(A{i},'DGBT') == 1 %Voltage Base
        j = 0;
        while true
            i = i + 1;
            if contains(A{i}(1:5),'99999') == 1
                break;
            end
            if strcmp(strtrim(A{i}(1)),'(') == 0 %
                j = j + 1;
                dgbt(j).grupo = A{i}(1:2);
                dgbt(j).tensao_base = str2double(A{i}(04:08));    %kV
            end
        end
    end
    if strfind(A{i},'FIM') == 1
        break;
    end
end

if isempty(dcte.base) == 1
    MVAbase = 100;
else
    MVAbase = dcte.base;
end
if isempty(dcte.dase) == 1
    MWbase = 100;
else
    MWbase = dcte.dase;
end

b = []; nMaxCar = 126;
for i = 1 : 1 : nMaxCar
    b = [b, ' '];
end

%Renumeração para o proj. com a prof. Elizandra da UFES
% for i = 1 : 1 : 107 
%     buses(i) = dbar(i).numero;
%     dbar(i).numero = i;    
% end
% cc = zeros(171,2);
% for i = 1 : 1 : 107
%     j = 108 - i;
%     for k = 1 : 1 : 171
%         if dlin(k).de == buses(j) && cc(k,1) == 0
%             dlin(k).de = j; cc(k,1) = 1;           
%         end
%         if dlin(k).para == buses(j) && cc(k,2) == 0
%             dlin(k).para = j;  cc(k,2) = 1;             
%         end        
%     end
% end

dataf = '0 /0 /0 ';

c = 0;
%%%%%%%%%%%Title Data%%%%%%%%%%
c = c + 1;
B{c} = b;
B{c}(2:9) = dataf; %Date, in format DD/MM/YY with leading zeros.
B{c}(11:30) = fullConvNum2Str(titu,11,30,c,'title');
t = num2str(MVAbase,'%.1f');
B{c}(32:37) = fullConvNum2Str(t,32,37,c,'MVABase');

%%%%%%%%%%%Bus Data%%%%%%%%%%
c = c+1;
nBus = length(dbar);
B{c} = b;
B{c}(1:16) = 'BUS DATA FOLLOWS';
nItens = string(strcat(num2str(nBus), {' '}, 'ITEMS'));
x2 = strlength(nItens);
B{c}(40:39+x2) = nItens;
for i = 1 : 1 : nBus
    c = c + 1;
    B{c} = b;
    %Columns  1- 4   Bus number (I) *
    t = num2str(dbar(i).numero);
    B{c}(1:4) = fullConvNum2Str(t,1,4,c,'bus number');
    %Columns  7-17   Name (A) (left justify) *
    t = num2str(dbar(i).nome);
    B{c}(7:17) = fullConvNum2Str(t,7,17,c,'name');
    %Columns 19-20   Load flow area number (I) Don't use zero! *
    t = num2str(dbar(i).area);
    B{c}(19:20) = fullConvNum2Str(t,19,20,c,'area number');
    %Columns 21-23   Loss zone number (I)
    %Columns 25-26   Type (I) *
    %                 0 - Unregulated (load, PQ)
    %                 1 - Hold MVAR generation within voltage limits, (PQ)
    %                 2 - Hold voltage within VAR limits (gen, PV)
    %                 3 - Hold voltage and angle (swing, V-Theta) (must always
    %                      have one)
    if isnan(dbar(i).tipo) == 1 || dbar(i).tipo == 0
        B{c}(26) = '0';
    elseif  dbar(i).tipo == 1
        B{c}(26) = '2';
    elseif  dbar(i).tipo == 3
        B{c}(26) = '1';
    elseif dbar(i).tipo == 2
        B{c}(26) = '3';
    end
    %Columns 28-33   Final voltage, p.u. (F) *
    t = num2str(dbar(i).tensao/1000,'%.3f');
    B{c}(28:33) = fullConvNum2Str(t,28,33,c,'final voltage');
    %Columns 34-40   Final angle, degrees (F) *
    if isnan(dbar(i).angulo) == 0
        t = num2str(dbar(i).angulo);
        B{c}(34:40) = fullConvNum2Str(t,34,40,c,'final angle');
    else
        B{c}(38:40) = '0.0';
    end
    %Columns 41-49   Load MW (F) *
    if isnan(dbar(i).demanda_ativa) == 0
        t = num2str(dbar(i).demanda_ativa);
        B{c}(41:49) = fullConvNum2Str(t,41,49,c,'load MW ');
    else
        B{c}(47:49) = '0.0';
    end
    %Columns 50-59   Load MVAR (F) *
    if isnan(dbar(i).demanda_reativa) == 0
        t = num2str(dbar(i).demanda_reativa);
        B{c}(50:59) = fullConvNum2Str(t,50,59,c,'load MVAR');
    else
        B{c}(57:59) = '0.0';
    end
    %Columns 60-67   Generation MW (F) *
    if isnan(dbar(i).potencia_ativa) == 0
        t = num2str(dbar(i).potencia_ativa);
        B{c}(60:67) = fullConvNum2Str(t,60,67,c,'generation MW');
    else
        B{c}(65:67) = '0.0';
    end
    %Columns 68-75   Generation MVAR (F) *
    if isnan(dbar(i).potencia_reativa) == 0
        t = num2str(dbar(i).potencia_reativa);
        B{c}(68:75) = fullConvNum2Str(t,68,75,c,'generation MVAR');
    else
        B{c}(73:75) = '0.0';
    end
    %Columns 77-83   Base KV (F)
    if dbar(i).grupo_base_tensao(1) ~= ' ' && dbar(i).grupo_base_tensao(2) ~= ' '
        n = length(dgbt);
        for j = 1 : 1 : n
            if strcmp(dbar(i).grupo_base_tensao,dgbt(j).grupo) == 1
                t = num2str(dgbt(j).tensao_base);
                B{c}(77:83) = fullConvNum2Str(t,77,83,c);
            end
        end
    else
        B{c}(81:83) = '1.0';
    end


    %Columns 85-90   Desired volts (pu) (F) (This is desired remote voltage if
    %                this bus is controlling another bus.
    %Columns 91-98   Maximum MVAR or voltage limit (F)
    if isnan(dbar(i).potencia_reativa_maxima) == 0
        t = num2str(dbar(i).potencia_reativa_maxima);
        B{c}(91:98) = fullConvNum2Str(t,91,98,c,'maximum MVAR');
    else
        B{c}(96:98) = '0.0';
    end
    %Columns 99-106  Minimum MVAR or voltage limit (F)
    if isnan(dbar(i).potencia_reativa_minima) == 0
        t = num2str(dbar(i).potencia_reativa_minima);
        B{c}(99:106) = fullConvNum2Str(t,99,106,c,'minimum MVAR ');
    else
        B{c}(104:106) = '0.0';
    end
    %Columns 107-114 Shunt conductance G (per unit) (F) *
    B{c}(112:114) = '0.0'; %Sem especificacao no ANAREDE
    %Columns 115-122 Shunt susceptance B (per unit) (F) *
    if isnan(dbar(i).shunt_barra) == 0
        t = num2str(dbar(i).shunt_barra/MVAbase);
        B{c}(115:122) = fullConvNum2Str(t,115,122,c,'shunt susceptance');
    else
        B{c}(120:122) = '0.0';
    end
    %Columns 124-127 Remote controlled bus number
    if isnan(dbar(i).barra_controlada) == 0
        t = num2str(dbar(i).shunt_barra/MVAbase);
        B{c}(124:127) = fullConvNum2Str(t,124,127,c,'controlled bus number');
    else
        B{c}(125:127) = '0.0';
    end
end
c = c+1;
B{c}(1:4) =  '-999';

%%%%%%%%%%%Branch data cards%%%%%%%%%%
c = c+1;
nBranches = length(dlin);
B{c} = b;
B{c}(1:19) = 'BRANCH DATA FOLLOWS';
nItens = string(strcat(num2str(nBranches), {' '}, 'ITEMS'));
x2 = strlength(nItens);
B{c}(40:39+x2) = nItens;
for i = 1 : 1 : nBranches
    c = c + 1;
    B{c} = b;
    %Columns  1- 4   Tap bus number (I) *
    t = num2str(dlin(i).de);
    B{c}(1:4) = fullConvNum2Str(t,1,4,c,'tap bus number');
    %Columns  6- 9   Z bus number (I) *
    t = num2str(dlin(i).para);
    B{c}(6:9) = fullConvNum2Str(t,6,9,c,'Z bus number');
    %Columns 11-12   Load flow area (I)
    %Columns 13-14   Loss zone (I)
    %Column  17      Circuit (I) * (Use 1 for single lines)
    %Column  19      Type (I) *
    %Columns 20-29   Branch resistance R, per unit (F) *
    if isnan(dlin(i).resistencia) == 0
        t = num2str(dlin(i).resistencia/100);
        B{c}(20:29) = fullConvNum2Str(t,20,29,c,'branch resistance');
    else
        B{c}(27:29) = '0.0';
    end
    %Columns 30-40   Branch reactance X, per unit (F) * No zero impedance lines
    if isnan(dlin(i).reatancia) == 0
        t = num2str(dlin(i).reatancia/100);
        B{c}(30:40) = fullConvNum2Str(t,30,40,c,'branch reactance');
    else
        B{c}(38:40) = '0.0';
    end
    %Columns 41-50   Line charging B, per unit (F) * (total line charging, +B)
    if isnan(dlin(i).susceptancia) == 0
        t = num2str(dlin(i).susceptancia/MVAbase);
        B{c}(41:50) = fullConvNum2Str(t,41,50,c,'line charging B');
    else
        B{c}(48:50) = '0.0';
    end
    %Columns 51-55   Line MVA rating No 1 (I) Left justify!
    %Columns 57-61   Line MVA rating No 2 (I) Left justify!
    %Columns 63-67   Line MVA rating No 3 (I) Left justify!
    %Columns 69-72   Control bus number
    %Column  74      Side (I)
    %                 0 - Controlled bus is one of the terminals
    %                 1 - Controlled bus is near the tap side
    %                 2 - Controlled bus is near the impedance side (Z bus)
    %Columns 77-82   Transformer final turns ratio (F)
    if isnan(dlin(i).tap) == 0
        t = num2str(dlin(i).tap);
        B{c}(77:82) = fullConvNum2Str(t,77,82,c,'transformer turns ratio');
    else
        B{c}(80:82) = '0.0';
    end
    %Columns 84-90   Transformer (phase shifter) final angle (F)
    if isnan(dlin(i).tap_defasagem) == 0
        t = num2str(dlin(i).tap_defasagem*pi/180);
        B{c}(84:90) = fullConvNum2Str(t,77,82,c,'transformer final angle');
    else
        B{c}(88:90) = '0.0';
    end
    %Columns 91-97   Minimum tap or phase shift (F)
    %Columns 98-104  Maximum tap or phase shift (F)
    %Columns 106-111 Step size (F)
    %Columns 113-119 Minimum voltage, MVAR or MW limit (F)
    %Columns 120-126 Maximum voltage, MVAR or MW limit (F)
end
c = c+1;
B{c}(1:4) =  '-999';

c = c+1;
B{c}(1:11) =  'END OF DATA';

nomeArquivo = strcat('.\Sistemas\',system,'.cdf');
%nomeArquivo = strcat('.\Sistemas\',system,'cdf.txt');
fout = fopen(nomeArquivo, 'wt+');

for i = 1:numel(B)
    if B{i} == -1
        fprintf(fout,'%s', B{i});
        break
    else
        fprintf(fout,'%s\n', B{i});
    end
end

fclose(fout);


function [resConv] = fullConvNum2Str(t,x1,x2,c, fname)
nCarc1 = length(t);
nCarc2 = x2-x1+1;
resConv = [];
if nCarc2 > nCarc1
    difCar = nCarc2 - nCarc1;
    for k = 1 : 1 : difCar
        resConv = [resConv, ' '];
    end
    for k = 1 : 1 : nCarc1
        resConv = [resConv, t(k)];
    end
elseif nCarc2 == nCarc1
    resConv = t;
else
    fprintf('Excesso de caracteres na linha %d do cartao - campo: %s \n', c, fname);
    for k = 1 : 1 : nCarc2
        resConv = [resConv, t(k)];
    end
end
end



