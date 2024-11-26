from numpy import zeros, ones, conj, exp, diag, asmatrix, asarray, concatenate, sum, linalg, sqrt, random, mean, var
from scipy.sparse import issparse, csr_matrix as sparse
from copy import deepcopy

import array as arr
import numpy as np
import math
import array


from calcYbus import calcYbus
from importSystemData import importSystemData
from runDirectMethod import runDirectMethod

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#Leitura do arquivo 
system = '3bus_modified.txt'
[dbus,dbranch,nbus,MVAbase] = importSystemData(system)

#Calculo da matriz Ybus
ybus = calcYbus(nbus, dbus, dbranch)    
sYbus = sparse(ybus)

#Dados do sistema
typeArr = arr.array('i',dbus['type'].values)
pl0 = dbus['lMW'].values / MVAbase
ql0  = dbus['lMVAR'].values / MVAbase
pg0 = dbus['gMW'].values / MVAbase
qg0  = dbus['gMVAR'].values / MVAbase
magV = dbus['magV'].values 
angV = dbus['angV'].values*math.pi/180

pv = [i for i in range(len(typeArr)) if typeArr[i] == 2]; npv = len(pv)
pq = [i for i in range(len(typeArr)) if typeArr[i] == 0]; npq = len(pq)
ref = [i for i in range(len(typeArr)) if typeArr[i] == 3]

nbus = npv + npq + 1

#Inicializacao das variaveis do método direto
wP = zeros((nbus,1)); wP[pv+pq] = 1 #autovalores de dP
wQ = zeros((nbus,1)); wQ[pq] = 1 #autovalores de dQ
t = zeros((1,1)) #tamanho do crescimento 
itMax = 20
tol = 1e-4

#################Caso base################# 
kp = ones((1,3))
[tload, dPtotal, it] = runDirectMethod(wP, wQ, tol, pq, npq, pv, npv, sYbus, magV, angV, pg0, qg0, pl0, ql0, nbus, itMax, kp)
tload_casobase = tload[0][0]
PdCB = (1+tload_casobase)*pl0
print('---------------RESULTADO DO CASO BASE---------------\n')
print('Acrescimo de carga total igual a ', dPtotal, 'pu', 'com lambda igual a ', tload_casobase)

#Sorteio aleatorio 
rng = np.random.default_rng(seed=None)

##################Algoritmo de entropia cruzada################# 
#Parametros do algoritmo 
nsamplesEC = 50
nElite = 10
PdElite = zeros((nElite,nbus))
itMaxEC = 100; sig2_min = 1e-3
t = 0.8; v = 0.2 #Fatores de suavizacao 
pos = np.where(pl0 == 0); pos = pos[0]; npos = np.prod(np.shape(pos)) #Vetor com os indices das barras sem cargas
pos2 = np.where(pl0 > 0); pos2 = pos2[0]; npos2 = np.prod(np.shape(pos2)) #Vetor com os indices das barras com cargas

#Parametros do sistema 
mu = ones((nbus,1))
sig2 = 0.025*ones((nbus,1))

#Exlcusao de direcao de crescimento para as barras sem carga
for i in range(npos):
    mu[pos[i],0] = 0
    sig2[pos[i],0] = 0

itEC = 0; percViol = 4/100; atingiu = 0; 
X = zeros((nbus,nsamplesEC))
fEC = zeros((nsamplesEC,3))
while True:
    itEC = itEC + 1    
    #Passo 1 - Criar o vetor de amostras X e avaliar a funcao objetivo    
    for i in range(nsamplesEC):     
        while True: 
            X[:,[i]] = mu[:,] + sqrt(sig2[:,]*abs(random.randn(nbus,1)))
            for j in range(npos):
                X[pos[j],i] = 0
            kp = X[:,[i]].reshape(nbus,)
            [tload, dPtotal, it] = runDirectMethod(wP, wQ, tol, pq, npq, pv, npv, sYbus, magV, angV, pg0, qg0, pl0, ql0, nbus, itMax, kp)
            fpen = 0
            if tload > 0 and it < itMax:
                fEC[i,0] = dPtotal + fpen
                fEC[i,1] = tload[0][0]
                fEC[i,2] = dPtotal
                break

    #Passo 2 - Criar o conjunto ordenado por S(X) que contém as amostras de elite 
    idxElite = np.argsort(fEC[:,0])
    Xel = deepcopy(X[:,idxElite[0:nElite]])
    
    #Passo 3 - Identificar se ha violacao da margem de seguranca 
    nViol = 0
    for i in range(nElite):
        if fEC[idxElite[i],1] < percViol:
            nViol += 1
    if atingiu == 0 and nViol >= 0.8*nElite:
        atingiu = 1
        print('\n---------------RESULTADO DA ENTROPIA CRUZADA---------------\n')
        print('\n Atingiu a condicao de violacao de ', percViol, 'de margem de carga em ', itEC, 'iteracoes \n')
        muViolEC = mu 
        sig2ViolEC = sig2

    #Passo 4 - Recaular a media e a variancia 
    mu0 = deepcopy(mu); sig20 = deepcopy(sig2)
    for i in range(nbus):        
        mu[i,0] = mean(Xel[i,:])
        sig2[i,0] = var(Xel[i,:])

    #Passo 5 - Aplicar fatores de suavizacao 
    mu = t*mu + (1-t)*mu0
    sig2 = v*sig2 + (1-v)*sig2

    #Passo 6 - Verificar condicao de parada 
    if itEC == itMaxEC or max(sig2) < sig2_min:
        if atingiu == 0:
            print('\n---------------RESULTADO DA ENTROPIA CRUZADA---------------\n')
        print('Acrescimo de carga total igual a ', fEC[idxElite[1],2], 'com lambda igual a ', fEC[idxElite[1],1],'\n')
        print('Valor de kp*lambda:\n')     
        for i in range(nbus):
            print('Barra', dbus['num'][i], ' ', fEC[idxElite[1],1]*X[i,idxElite[1]])
        for i in range(nElite):
            for j in range(nbus):
                PdElite[i,j] = (1 + fEC[idxElite[i],1]*Xel[j,i])*pl0[j]
        break


PdEC = zeros((nbus,1))
for j in range(nbus):
    PdEC[j,0] = (1 + fEC[idxElite[0],1]*Xel[j,0])*pl0[j]
    

##################Algoritmo de evolucao diferencial################# 
#Parametros iniciais
itMaxED = 100 #Numero maximo de iteracoes
LP = 20 #Quantidade de iteracoes no periodo de aprendizagem
nEstrategias = 4 #Quantidade de mutacoes avaliadas
cont = array.array('i',(0 for i in range(0,nEstrategias))) 
CRm = 0.5*ones((nEstrategias,1)) #Taxa de crossover inicial
p = 1/nEstrategias*ones((1,nEstrategias)) #Inicializa as probabilidades de escolher cada mutacao
ns = array.array('i',(0 for i in range(0,nEstrategias))) #Numero de mutacoes bem sucedidas por estrategia 
nf = array.array('i',(0 for i in range(0,nEstrategias))) #Numero de falhas de cada uma das estrategias
txns = zeros((1,nEstrategias))
prob = zeros((1,nEstrategias))
#Populacao inicial (resultado da entropia cruzada)
pop = deepcopy(X[:,idxElite[0:nElite]])
nInd = nElite
fpop = deepcopy(fEC[0:nElite,:])
estrategia = array.array('i',(0 for i in range(0,nInd)))

CRMemory = zeros((nInd,1))
CR = zeros((1,nEstrategias))
sel = array.array('i',(0 for i in range(0,5)))
itED = 0
while True:
    itED += 1
    if itED == 1:
        novaPop = deepcopy(pop)
        fnovaPop = deepcopy(fpop)
    posMelhorInd = np.argmin(fpop[:,0])
    Z = zeros((nInd,1))
    for i in range(nInd):
        #Fator de mutacao
        F = random.normal(0.5, 0.1) #Sorteio da gaussiana com mu = 0,5 e sig = 0.1
        #Selecao dos individuos 
        aux = random.permutation(nElite)
        j = 0
        k = 0
        while True:
            if i != aux[j]:              
                sel[k] = aux[j]
                k += 1
            j += 1
            if k == 5:
                break
        #Selecao da estregia de mutacao durante o periodo de aprendizagem 
        if itED <= LP:
            estrategia[i] = random.randint(nEstrategias)
        else:
            estrategia[i] = estrategiaFinal    
        jrand = random.randint(npos2)
        if estrategia[i] == 0:
            while True:
                pc = random.normal(CRm[0,0], 0.1)
                if pc > 0 and pc < 1:
                    break
            CRMemory[i] = pc
            for j in range(npos2):               
                if random.rand() < pc or j == jrand:
                    novaPop[pos2[j],i] = pop[pos2[j],sel[0]]+ F*(pop[pos2[j],sel[1]] - pop[pos2[j],sel[2]])
                else:
                    novaPop[pos2[j],i] = pop[pos2[j],i]
        elif estrategia[i] == 1:
            while True:
                pc = random.normal(CRm[1,0], 0.1)
                if pc > 0 and pc < 1:
                    break    
            CRMemory[i] = pc
            for j in range(npos2):               
                if random.rand() < pc or j == jrand:
                    novaPop[pos2[j],i] = pop[pos2[j],i] + F*(pop[pos2[j],posMelhorInd]-pop[pos2[j],i]) + F*(pop[pos2[j],sel[0]]-pop[pos2[j],sel[1]]) + F*(pop[pos2[j],sel[2]]-pop[pos2[j],sel[3]])
                else:
                    novaPop[pos2[j],i] = pop[pos2[j],i]
        elif estrategia[i] == 2:
            while True:
                pc = random.normal(CRm[2,0], 0.1)
                if pc > 0 and pc < 1:
                    break    
            CRMemory[i] = pc
            for j in range(npos2): 
                if random.rand() < pc or j == jrand:
                    novaPop[pos2[j],i] = pop[pos2[j],i] + F*(pop[pos2[j],sel[1]]-pop[pos2[j],sel[2]]) + F*(pop[pos2[j],sel[3]]-pop[pos2[j],sel[4]]) 
                else:
                    novaPop[pos2[j],i] = pop[pos2[j],i]
        elif estrategia[i] == 3:
            F2 = random.rand()
            for j in range(npos2):
                novaPop[pos2[j],i] = pop[pos2[j],i] + F2*(pop[pos2[j],sel[0]]-pop[pos2[j],i]) + F*(pop[pos2[j],sel[1]]-pop[pos2[j],sel[2]])
        kp = novaPop[:,i]
        [tload, dPtotal, it] = runDirectMethod(wP, wQ, tol, pq, npq, pv, npv, sYbus, magV, angV, pg0, qg0, pl0, ql0, nbus, itMax, kp)
        fpen = 0
        fnovaPop[i,0] = dPtotal + fpen
        fnovaPop[i,1] = tload[0][0]
        fnovaPop[i,2] = dPtotal
        Z[i,0] = fnovaPop[i,0]
        if tload < 0 or it >= itMax or min(novaPop[:,i]) < 0:
            novaPop[:,i] = pop[:,i]
            fnovaPop[i,:] = fpop[i,:]
        if i == 9:
            teste = 2
    #Processo de selecao 
    for i in range(nInd):
        if fnovaPop[i,0] < fpop[i,0]:
            pop[:,i] = novaPop[:,i]
            fpop[i,:] = fnovaPop[i,:]  
            if itED <= LP:
                ns[estrategia[i]] += 1
                cont[estrategia[i]] += 1
            if estrategia[i] != 4: 
                aux = zeros((1,nEstrategias))
                aux[0,estrategia[i]] = CRMemory[i][0]
                CR = concatenate((CR,aux),axis = 1)
        else: 
            if itED <= LP:
                nf[estrategia[i]] += 1
    if itED == LP: 
        sumts = 0 #Somatorio das taxas de sucesso
        estrategiaFinal = 0
        for i in range(nEstrategias):
            txns[0][i] = ns[i] / (ns[i]+nf[i]) #Calculo da taxa de sucesso 
            sumts += txns[0][i]
        for i in range(nEstrategias):
            prob[0][i] = txns[0][i] / sumts
            if prob[0][i] > prob[0][estrategiaFinal]:
                estrategiaFinal = i
    if itED == itMaxED:
        posMelhorInd = np.argmin(fpop[:,0])
        print('\n---------------RESULTADO DO MÉTODO SADE---------------\n')
        print('Acrescimo de carga total igual a ', fpop[posMelhorInd,2], 'com lambda igual a ', fpop[posMelhorInd,1],'\n')
        print('Valor de kp*lambda:\n')     
        for i in range(nbus):
            print('Barra', dbus['num'][i], ' ', fpop[posMelhorInd][1]*pop[i,posMelhorInd])
        PdED = (1+fpop[posMelhorInd,1]*pop[:,posMelhorInd])*pl0
        break

            
teste = 1



import matplotlib.pyplot as plt
import scipy.io as scio

x = scio.loadmat('res_x.mat')
x = x['x']
y = scio.loadmat('res_y.mat')
y = y['y']

plt.figure(1)
plt.scatter(x,y, s=80, facecolors='none', edgecolors='b')

# %Caso base
pl0 = pl0.reshape(nbus,1)
PdCB = PdCB.reshape(nbus,1)
PdEC = PdEC.reshape(nbus,1)
PdED = PdED.reshape(nbus,1)
plt.scatter(pl0[1],pl0[2], s=80, facecolors='none', edgecolors='b')
plt.plot([pl0[1],PdCB[1]],[pl0[2],PdCB[2]],'k--')
plt.scatter(PdCB[1],PdCB[2], s=80, facecolors='none', edgecolors='k')
plt.plot([pl0[1][0],PdEC[1][0]],[pl0[2][0],PdEC[2][0]],'r')
plt.scatter(PdEC[1],PdEC[2], s=80, facecolors='none', edgecolors='r')
plt.plot([pl0[1][0],PdED[1][0]],[pl0[2],PdED[2]],'m')
plt.scatter(PdED[1],PdED[2], s=80, facecolors='none', edgecolors='m')

# xlim([0 1.2])
# ylim([0 2.2])

for i in range(nElite):
    plt.plot([pl0[1][0],PdElite[i][1]],[pl0[2][0],PdElite[i][2]],'r--')
    plt.scatter(PdElite[i][1],PdElite[i][2], s=80, facecolors='none', edgecolors='r')

plt.text(PdCB[1][0],PdCB[2][0],'On-scenario')
plt.text(PdEC[1][0],PdEC[2][0],'Minimum-\lambda of cross-entropy algorithm')
plt.text(PdED[1][0],PdED[2][0],'Minimum-\lambda of SADE algorithm')

plt.xlabel ('P_2(pu)')
plt.ylabel ('P_3(pu)')

plt.grid('minor')

plt.show()
# % set(gca,'FontSize',15,'FontName','Times')
# % set(gcf,'Paperunits','inches','PaperPosition',[0 0 20 15],'PaperSize',[20 15])
# % pause(1);
# % print('Caso2','-dpdf','-r400')


   



