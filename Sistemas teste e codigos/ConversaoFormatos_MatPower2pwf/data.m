function [NB,NA,A,pq,pv,ref,gen,Ybus,BC,V,teta,Pg,Pmax,Qg,Pd,Qd,Qmax,Qmin, ...
            Vmax,Vmin,Sb,names,buses,branches] = data(casobase)

%%DATA_FULL Import system data and build Ybus
if nargin==0
    casobase = case3; %O default é o sistema IEEE 14bus
end

mpc = casobase; %Importação do caso escolhido


%% Importar os ramos
branches = mpc.branch;

%% Importar os nomes das barras
names = mpc.bus_name;
BC = mpc.bus(:,6);

%% Definição dos parâmetros dimensionais do case
NB = length(mpc.bus(:,1)); %Número de barras
NL = length(mpc.branch(:,1)); %Número de ramos
%% Dados das áreas
A = mpc.bus(:,7); NA = max(mpc.bus(:,7)); %Número de áreas
buses = mpc.bus(:,1); %Barras originais
% pqo = mpc.bus((mpc.bus(:,2)==1),1); %Barras PQ
% pvo = mpc.bus((mpc.bus(:,2)==2),1); %Barras PV
% refo = mpc.bus((mpc.bus(:,2)==3),1); %Barra de referência
% busesgen = sort([refo;pvo]);
pq = find(mpc.bus(:,2)==1); %Barras PQ 
pv = find(mpc.bus(:,2)==2); %Barras PV
ref = find(mpc.bus(:,2)==3); %Barra de referência
Sb = mpc.baseMVA; %Potência base do sistema
gen = sort([ref;pv]); %Barras com geração

%% Limites de reativo das barras de geração, já em p.u:
Qmax = mpc.gen(:,4)/Sb; Qmin = mpc.gen(:,5)/Sb;
%Qmax = Qmax(Qmax<5000); Qmin = Qmin(Qmin>-5000);

%% Limites de tensão de cada uma das barras, já em p.u:
Vmax = mpc.bus(:,12); Vmin = mpc.bus(:,13);

%% Gerações das barras 
Pmax = zeros(NB,1); Pmax(gen) = mpc.gen(:,9);
Pg = zeros(NB,1); Pg(gen) = mpc.gen(:,2);
Qg = zeros(NB,1); Qg(gen) = mpc.gen(:,3);

%% Demandas e tensões fasoriais das barras      
V = mpc.bus(:,8); 
V(gen) = mpc.gen(:,6);
teta = mpc.bus(:,9); 
Pd = mpc.bus(:,3); Qd = mpc.bus(:,4);

%% Transformar todas as potências para p.u
esc = 1;
Pg = esc*Pg/Sb; Qg = Qg/Sb;
Pmax = Pmax/Sb;
Pd = esc*Pd/Sb; Qd = esc*Qd/Sb;

%% Transformar os ângulos para radianos
teta = teta*pi/180;

%% Computar a matriz de admitância
% Dados dos ramos modificados
r = mpc.branch(:,3); x = mpc.branch(:,4); b = mpc.branch(:,5);
status = mpc.branch(:,11);
tap = ones(NL, 1);              %% default tap ratio = 1
i = find(mpc.branch(:,9));      %% indices of non-zero tap ratios
tap(i) = mpc.branch(i,9);       %% assign non-zero tap ratios
tap = tap .* exp(1j*pi/180 * mpc.branch(:,10)); %% add phase shifters

%% Computar os elementos série e shunt
Ys = status ./ (r + 1i * x);  %% Adimitância série
Bc = status .* b;            %% Elemento shunt da barra
Ysh = (mpc.bus(:,5)+1i*mpc.bus(:,6))/Sb;

%% Montagem da matriz
Ytt = Ys + 1i*Bc/2;
Yff = Ytt ./ (tap .* conj(tap));
Yft = - Ys ./ conj(tap);
Ytf = - Ys ./ tap;

%% Indices dos ramos
f = mpc.branch(:,1);  % "from" buses
t = mpc.branch(:,2);    % "to" buses
fn = zeros(size(f)); tn = zeros(size(t));
for j=1:NL
    fj = mpc.branch(j,1); tj = mpc.branch(j,2);
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