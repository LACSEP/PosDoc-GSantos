clc; clear;

addpath('.\Resultados\');

fMC1 = load('SMC_rho0.0_3e3.mat'); 
fMC1 = fMC1.S_MC; fMC1 = fMC1(:,1);
j = 0;
for i = 1 : 1 : 3000
    if fMC1(i) < 0.07
        j = j + 1;        
    end
    p(i) = j/i;
end
pMCS(1,1) = p(i)*100;
eMCS(1,1) = 2.5758*std(p)/sqrt(3000)/mean(p)*100; 

fMC2 = load('SMC_rho0.4_3e3.mat'); 
fMC2 = fMC2.S_MC; fMC2 = fMC2(:,1);
j = 0;
for i = 1 : 1 : 3000
    if fMC2(i) < 0.07
        j = j + 1;        
    end
    p(i) = j/i;
end
pMCS(2,1) = p(i)*100;
eMCS(2,1) = 2.5758*std(p)/sqrt(3000)/mean(p)*100; 

fMC3 = load('SMC_rho0.8_3e3.mat'); 
fMC3 = fMC3.S_MC; fMC3 = fMC3(:,1);
j = 0;
for i = 1 : 1 : 3000
    if fMC3(i) < 0.07
        j = j + 1;        
    end
    p(i) = j/i;
end
pMCS(3,1) = p(i)*100;
eMCS(3,1) = 2.5758*std(p)/sqrt(3000)/mean(p)*100; 

fMC4 = load('SMC_rho0.0_1e4.mat'); 
fMC4 = fMC4.S_MC; fMC4 = fMC4(:,1);
j = 0;
for i = 1 : 1 : 10000
    if fMC4(i) < 0.07
        j = j + 1;        
    end
    p(i) = j/i;
end
pMCS(1,2) = p(i)*100;
eMCS(1,2) = 2.5758*std(p)/sqrt(10000)/mean(p)*100; 

fMC5 = load('SMC_rho0.4_1e4.mat'); 
fMC5 = fMC5.S_MC; fMC5 = fMC5(:,1);
j = 0;
for i = 1 : 1 : 10000
    if fMC5(i) < 0.07
        j = j + 1;        
    end
    p(i) = j/i;
end
pMCS(2,2) = p(i)*100;
eMCS(2,2) = 2.5758*std(p)/sqrt(10000)/mean(p)*100; 

fMC6 = load('SMC_rho0.8_1e4.mat'); 
fMC6 = fMC6.S_MC; fMC6 = fMC6(:,1);
j = 0;
for i = 1 : 1 : 10000
    if fMC6(i) < 0.07
        j = j + 1;        
    end
    p(i) = j/i;
end
pMCS(3,2) = p(i)*100;
eMCS(3,2) = 2.5758*std(p)/sqrt(10000)/mean(p)*100; 

fMC7 = load('SMC_rho0.0.mat'); 
fMC7 = fMC7.S_MC; fMC7 = fMC7(:,1);
j = 0;
for i = 1 : 1 : 1e6
    if fMC7(i) < 0.07
        j = j + 1;        
    end
    p(i) = j/i;
end
pMCS(1,3) = p(i)*100;
eMCS(1,3) = 2.5758*std(p)/sqrt(1e6)/mean(p)*100; 

fMC8 = load('SMC_rho0.4.mat'); 
fMC8 = fMC8.S_MC; fMC8 = fMC8(:,1);
j = 0;
for i = 1 : 1 : 1e6
    if fMC8(i) < 0.07
        j = j + 1;        
    end
    p(i) = j/i;
end
pMCS(2,3) = p(i)*100;
eMCS(2,3) = 2.5758*std(p)/sqrt(1e6)/mean(p)*100; 

fMC9 = load('SMC_rho0.8.mat'); 
fMC9 = fMC9.S_MC; fMC9 = fMC9(:,1);
j = 0;
for i = 1 : 1 : 1e6
    if fMC9(i) < 0.07
        j = j + 1;        
    end
    p(i) = j/i;
end
pMCS(3,3) = p(i)*100;
eMCS(3,3) = 2.5758*std(p)/sqrt(1e6)/mean(p)*100; 

p = [];
fCE1 = load('SEC_rho0.0.mat');
fCE1 = fCE1.S2_EC; fCE1 = fCE1(:,1);
wCE1 = load('WEC_rho0.0.mat');
wCE1 = wCE1.W;
j = 0;
for i = 1 : 1 : 3000
    if fCE1(i) < 0.07
        I1(i) = 1;
    else
        I1(i) = 0;
    end
    p(i) = sum(I1.*wCE1(1:i))/i;    
end
pCE(1,1) = p(i)*100;
eCE(1,1) = 2.5758*std(p)/sqrt(3000)/mean(p)*100; 

p = []; I1 = [];
fCE2 = load('SEC_rho0.4.mat');
fCE2 = fCE2.S2_EC; fCE2 = fCE2(:,1);
wCE2 = load('WEC_rho0.4.mat');
wCE2 = wCE2.W;
j = 0;
for i = 1 : 1 : 3000
    if fCE2(i) < 0.07
        I1(i) = 1;
    else
        I1(i) = 0;
    end
    p(i) = sum(I1.*wCE2(1:i))/i;    
end
pCE(2,1) = p(i)*100;
eCE(2,1) = 2.5758*std(p)/sqrt(3000)/mean(p)*100; 

p = []; I1 = [];
fCE3 = load('SEC_rho0.8.mat');
fCE3 = fCE3.S2_EC; fCE3 = fCE3(:,1);
wCE3 = load('WEC_rho0.8.mat');
wCE3 = wCE3.W;
j = 0;
for i = 1 : 1 : 3000
    if fCE3(i) < 0.07
        I1(i) = 1;
    else
        I1(i) = 0;
    end
    p(i) = sum(I1.*wCE3(1:i))/i;    
end
pCE(3,1) = p(i)*100;
eCE(3,1) = 2.5758*std(p)/sqrt(3000)/mean(p)*100; 