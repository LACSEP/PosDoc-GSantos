%Formato dos arquivos: https://labs.ece.uw.edu/pstca/formats/cdf.txt
%Dados dos sistemas: https://labs.ece.uw.edu/pstca/

function [NB,NA,areas,pq,pv,ref,gen,Ybus,V,teta,Pg,Qg,Pd,Qd,Qmax,Qmin, ...
    Vmax,Vmin,Sb,names,buses] = systemData(system)

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

Sb = str2double(A{1}(32:37));

n = length(A);
i = 0;
while true
    i = i + 1;
    if strfind(A{i},'BUS DATA FOLLOWS') == 1
        j = 0;
        while true
            i = i + 1;
            j = j + 1;
            if contains(A{i}(1:2),'-9') == 1 || contains(A{i}(1:3),'-99') == 1 || contains(A{i}(1:4),'-999') == 1
                break;
            end
            dbus(j).num = str2double(A{i}(1:4));
            dbus(j).name = A{i}(6:17);
            dbus(j).area = str2double(A{i}(20));
            dbus(j).type = str2double(A{i}(26));
            dbus(j).magV = str2double(A{i}(28:33));
            dbus(j).angV = str2double(A{i}(34:40));
            dbus(j).lMW = str2double(A{i}(42:49));
            dbus(j).lMVAR = str2double(A{i}(51:59));
            dbus(j).gMW = str2double(A{i}(61:67));
            dbus(j).gMVAR = str2double(A{i}(69:75));
            dbus(j).basekV = str2double(A{i}(78:83));
            dbus(j).maxMVAR = str2double(A{i}(92:98));
            dbus(j).minMVAR = str2double(A{i}(100:106));
            dbus(j).shuntG = str2double(A{i}(107:114));
            dbus(j).shuntB = str2double(A{i}(115:122));
        end
    end
    if strfind(A{i},'BRANCH DATA') == 1
        j = 0;
        while true
            i = i + 1;
            j = j + 1;
            if contains(A{i}(1:2),'-9') == 1 || contains(A{i}(1:3),'-99') == 1 || contains(A{i}(1:4),'-999') == 1
                break;
            end
            dbranch(j).iBus = str2double(A{i}(1:4));
            dbranch(j).jBus = str2double(A{i}(6:9));
            dbranch(j).area = str2double(A{i}(12));
            dbranch(j).type = str2double(A{i}(19));
            dbranch(j).lineR = str2double(A{i}(21:29));
            dbranch(j).lineX = str2double(A{i}(31:40));
            dbranch(j).lineB = str2double(A{i}(42:50));
            dbranch(j).magTap = str2double(A{i}(77:82));
            dbranch(j).angTap = str2double(A{i}(85:90));
        end
    end
    if strfind(A{i},'END OF DATA') == 1
        break;
    end
end
names = vertcat(dbus.name);
buses = vertcat(dbus.num);
NB = length(dbus);  %Número de barras
NL = length(dbranch); %%Número de ramos
areas = vertcat(dbus.area);
NA = max(areas); %Número de áreas
pq = find(vertcat(dbus.type)==0); %Barras PQ
pv = find(vertcat(dbus.type)==2); %Barras PV
ref = find(vertcat(dbus.type)==3); %Barra ref
gen = sort([ref;pv]); %Barras com geração

%% Limites de reativo das barras de geração, já em p.u:
Qmax = vertcat(dbus.maxMVAR)/Sb; Qmin = vertcat(dbus.minMVAR)/Sb;
Qmax = Qmax(Qmax<5000); Qmin = Qmin(Qmin>-5000);


%% Limites de tensão de cada uma das barras, já em p.u: (NAO INFORMADO NO CARTÃO DO IEEE!)
Vmax = zeros(NB,1);
Vmin = zeros(NB,1);

%% Gerações das barras (em p.u.)
Pg = vertcat(dbus.gMW)/Sb;
Qg = vertcat(dbus.gMVAR)/Sb;

%% Demandas e tensões fasoriais das barras (em p.u.)
V = vertcat(dbus.magV);
teta = vertcat(dbus.angV)*pi/180; %em radianos
Pd = vertcat(dbus.lMW)/Sb;
Qd = vertcat(dbus.lMVAR)/Sb;

%% Computar a matriz de admitância
r = vertcat(dbranch.lineR); x = vertcat(dbranch.lineX); b = vertcat(dbranch.lineB);
status = ones(NL,1); %(NAO INFORMADO NO CARTAO DO IEEE
tap = ones(NL, 1);              %% default tap ratio = 1
tapidx = find(vertcat(dbranch.magTap));      %% indices of non-zero tap ratios
aux = vertcat(dbranch.magTap);
tap(tapidx) = aux(tapidx); tap = tap .* exp(1i*pi/180 * vertcat(dbranch.angTap)); %% add phase shifters


%% Computar os elementos série e shunt
Ys = status ./ (r + 1i * x);  %% Adimitância série
Bc = status .* b;            %% Elemento shunt da barra
Ysh = vertcat(dbus.shuntG) +1i*vertcat(dbus.shuntB);

%% Montagem da matriz
Ytt = Ys + 1i*Bc/2;
Yff = Ytt ./ (tap .* conj(tap));
Yft = - Ys ./ conj(tap);
Ytf = - Ys ./ tap;

%% Indices dos ramos
f = vertcat(dbranch.iBus);  % "from" buses
t = vertcat(dbranch.jBus);    % "to" buses
fn = zeros(size(f)); tn = zeros(size(t));
for j=1:NL
    fj = f(j); tj = t(j);
    fn(j,1) = find(buses==fj); tn(j,1) = find(buses==tj);
end


%% Matrizes de conexão
Cf = sparse(1:NL, fn, ones(NL, 1), NL, NB);
Ct = sparse(1:NL, tn, ones(NL, 1), NL, NB);

%% Construir Yf e Yt
Yf = sparse(1:NL, 1:NL, Yff, NL, NL) * Cf + sparse(1:NL, 1:NL, Yft, NL, NL) * Ct;
Yt = sparse(1:NL, 1:NL, Ytf, NL, NL) * Cf + sparse(1:NL, 1:NL, Ytt, NL, NL) * Ct;

%% Construir Ybus
Ybus = Cf' * Yf + Ct' * Yt + sparse(1:NB, 1:NB, Ysh, NB, NB);
%Ybus = full(Ybus); %Converter pra matrix completa

end



