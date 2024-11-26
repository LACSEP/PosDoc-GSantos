clc; clear;

addpath('.\Sistemas')

casobase = case39;
filename = 'case39'; %nome do arquivo final .pwf


[nbus,NA,area,pq,pv,ref,gen,Ybus,BC,Vm,Va,Pg0,PMax,Qg,Pd0,Qd0,Qmax,Qmin,Vmax,Vmin,Sb,names,buses,branches] = data(casobase);

fpart = [];

% Fator de participacao para reserva de potência das maquinas sincronas que
% atuam como geradores sincronos (Pg0 ~= 0)
nger = length(gen);
fpart = zeros(nbus,1);
fpart_num = zeros(nger,1);
fpart_den = 0;
for i = 1 : 1 : nger
    if Pg0(gen(i)) ~= 0
        fpart_num(gen(i)) = PMax(gen(i))-Pg0(gen(i));
        fpart_den = fpart_den + fpart_num(gen(i));
    end
end
for i = 1 : 1 : nger
    if Pg0(gen(i)) ~= 0
        fpart(gen(i)) = fpart_num(gen(i))/fpart_den;       
    end
end

dirfile = strcat('.\Sistemas\',filename,'.pwf');

TITU = strcat('Sistema de',32, num2str(nbus),' barras - ', datestr(datetime));

%tipo de barra
type = zeros(nbus,1);
type(pv) = 1;
type(ref) = 2;

%potencia reativa minima / maxima
qmin = zeros(nbus,1); qmin(gen) = Qmin*Sb;
qmax = zeros(nbus,1); qmax(gen) = Qmax*Sb;

%circuitos de ramos 
nbranches = length(branches);
circuits = zeros(nbranches,1);
for i = 1 : 1 : nbranches
    if i == 1 
        circuits(i) = 1;
    else 
        if circuits(i) == 0
            circuits(i) = 1;            
            k = 1;
            j = i;
            while true                                 
                j = j + 1;
                if j > nbranches
                    break;
                end
                if branches(i,1) == branches(j,1) && branches(i,2) == branches(j,2)
                    k = k + 1;
                    circuits(j) = k;                    
                end
                if branches(i,1) == branches(j,2) && branches(i,2) == branches(j,1)
                    k = k + 1;
                    circuits(j) = k;                    
                end              
            end
        end
    end
end

b = []; nMaxCar = 120;
for i = 1 : 1 : nMaxCar
    b = [b, ' '];
end

c = 0;
%%%%%%%%%%%Title Data%%%%%%%%%%
c = c + 1;
B{c} = 'TITU';
c = c + 1;
B{c} = TITU;

%%%%%%%%%%%Codigo de Execução DBAR%%%%%%%%%%
c = c+1;
B{c} = b;
B{c}(1:4) = 'DBAR';
c = c+1;
B{c} = '(Num)OETGb(   nome   )Gl( V)( A)( Pg)( Qg)( Qn)( Qm)(Bc  )( Pl)( Ql)( Sh)Are(Vf)M(1)(2)(3)(4)(5)(6)(7)(8)(9)(10';
for i = 1 : 1 : nbus
    c = c+1;
    B{c} = b;
    t = num2str(buses(i));
    B{c}(1:5) = fullConvNum2Str(t,1,5,c,'numero da barra');
    t = num2str(type(i));
    B{c}(8) = fullConvNum2Str(t,8,8,c,'tipo da barra');
%     t = names{i};
%     B{c}(11:22) = fullConvNum2Str(t,11,22,c,'nome da barra');
    t = num2str(round(Vm(i)*1000));
    B{c}(25:28) = fullConvNum2Str(t,25,28,c,'magnitude da tensão');
    t = num2str(round(Va(i)*180/pi,1));
    B{c}(29:32) = fullConvNum2Str(t,29,32,c,'ângulo da tensão');
    t = num2str(round(Pg0(i)*Sb));
    B{c}(33:37) = fullConvNum2Str(t,33,37,c,'potência ativa gerada');
    t = num2str(round(Qg(i)*Sb));
    B{c}(38:42) = fullConvNum2Str(t,38,42,c,'potência reativa gerada');
    t = num2str(round(qmin(i)));
    B{c}(43:47) = fullConvNum2Str(t,43,47,c,'potência reativa mínima gerada');
    t = num2str(round(qmax(i)));
    B{c}(48:52) = fullConvNum2Str(t,48,52,c,'potência reativa mínima gerada');
    t = num2str(round(Pd0(i)*Sb,1));
    B{c}(59:63) = fullConvNum2Str(t,59,63,c,'potência ativa demandada');
    t = num2str(round(Qd0(i)*Sb,1));
    B{c}(64:68) = fullConvNum2Str(t,64,68,c,'potência reativa demandada');
    t = num2str(round(BC(i),1));
    B{c}(69:73) = fullConvNum2Str(t,69,73,c,'potência reativa demandada');
    t = num2str(area(i));
    B{c}(74:76) = fullConvNum2Str(t,74,76,c,'área da barra');
end
c = c+1;
B{c}(1:5) =  '99999';

%%%%%%%%%%%Codigo de Execução DLIN%%%%%%%%%%
c = c+1;
B{c} = 'DLIN';                                                                                                         
c = c+1;
B{c} = '(De )d O d(Pa )NcEPM( R% )( X% )(Mvar)(Tap)(Tmn)(Tmx)(Phs)(Bc  )(Cn)(Ce)Ns(Cq)(1)(2)(3)(4)(5)(6)(7)(8)(9)(10';
for i = 1 : 1 : nbranches
    c = c+1;
    B{c} = b;
    t = num2str(branches(i,1));
    B{c}(1:5) = fullConvNum2Str(t,1,5,c,'número barra de');
    t = num2str(branches(i,2));
    B{c}(11:15) = fullConvNum2Str(t,1,5,c,'número barra para');
    t = num2str(circuits(i));
    B{c}(16:17) = fullConvNum2Str(t,16,17,c,'número circuito');
    t = num2str(round(branches(i,3)*100,3));
    B{c}(21:26) = fullConvNum2Str(t,21,26,c,'resistência da linha');
    t = num2str(round(branches(i,4)*100,3));
    B{c}(27:32) = fullConvNum2Str(t,27,32,c,'reatância da linha');
    t = num2str(round(branches(i,5)*100,3));
    B{c}(33:38) = fullConvNum2Str(t,27,32,c,'susceptância da linha');
    if branches(i,9)~= 0
        t = num2str(round(branches(i,9),3));
        B{c}(39:43) = fullConvNum2Str(t,39,43,c,'tap do transformador');
        t = num2str(round(branches(i,10),3));
        B{c}(54:58) = fullConvNum2Str(t,39,43,c,'defasagem do transformador');
    end
    if branches(i,6)~= 0
        t = num2str(branches(i,6));
        B{c}(65:68) = fullConvNum2Str(t,65,68,c,'capacidade normal da linha');
    end
    if branches(i,8)~= 0
        t = num2str(branches(i,8));
        B{c}(69:72) = fullConvNum2Str(t,69,72,c,'capacidade em emergência da linha');
    end
end
c = c+1;
B{c}(1:5) =  '99999';

if isempty(fpart) == 0 %Fator de participacao foi informado
    %%%%%%%%%%%Codigo de Execução DGER%%%%%%%%%%
    c = c+1;
    B{c} = 'DGER';   
    c = c+1;
    B{c} = 'DGER';   
    B{c} = '(No ) O (Pmn ) (Pmx ) ( Fp) (FpR) (FPn) (Fa) (Fr) (Ag) ( Xq) (Sno) (Est)';
    for i = 1 : 1 : nger
        c = c+1;
        j = gen(i);
        B{c} = b;
        t = num2str(buses(j));
        B{c}(1:5) = fullConvNum2Str(t,1,5,c,'numero do gerador');
        t = num2str(round(fpart(j)*100,2));
        B{c}(23:27) = fullConvNum2Str(t,23,27,c,'fator de participacao');
    end
end
c = c+1;
B{c}(1:5) =  '99999';

c = c+1;
B{c}(1:3) =  'FIM';

fout = fopen(dirfile, 'wt+');

for i = 1:numel(B)    
    if B{i} == -1
        fprintf(fout,'%s', B{i});
        break
    else
        fprintf(fout,'%s\n', B{i});
    end
end


fclose(fout);




