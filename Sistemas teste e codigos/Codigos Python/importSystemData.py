import os
from os.path import dirname, exists, realpath
from pandas import DataFrame as DF


def importSystemData(system):
    #Inicializacao da leitura
    systemDir = os.path.join(os.path.dirname(__file__), 'Sistemas', system)    
    f = open(f'{systemDir}', 'r', encoding='latin-1')
    lines = f.readlines()
    f.close()
    pwf2py = {}

    MVAbase = float(lines[0][31:37].strip())

    dbus = dict()
    dbus['num'] = list()
    dbus['name'] = list()
    dbus['area'] = list()
    dbus['type'] = list()
    dbus['magV'] = list()
    dbus['angV'] = list()
    dbus['lMW'] = list()
    dbus['lMVAR'] = list()
    dbus['gMW'] = list()
    dbus['gMVAR'] = list()
    dbus['basekV'] = list()
    dbus['maxMVAR'] = list()
    dbus['minMVAR'] = list()
    dbus['shuntG'] = list()
    dbus['shuntB'] = list()

    dbranch = dict()
    dbranch['iBus'] = list()
    dbranch['jBus'] = list()
    dbranch['area'] = list()
    dbranch['type'] = list()
    dbranch['lineR'] = list()
    dbranch['lineX'] = list()
    dbranch['lineB'] = list()
    dbranch['magTap'] = list()
    dbranch['angTap'] = list()

    endBlock = ('-9','-99', '-999')

    i = 0 
    while lines[i].strip() != 'END OF DATA': #Loop de leitura dos arquivos .txt 
        i += 1 
        if 'BUS DATA FOLLOWS' in lines[i].strip():        
            i += 1
            nbus = 0
            while lines[i][0:4].strip() not in endBlock: 
                dbus['num'].append(lines[i][0:4].strip())
                dbus['name'].append(lines[i][5:17].strip())
                dbus['area'].append(lines[i][19:20].strip())
                dbus['type'].append(lines[i][25:26].strip())
                dbus['magV'].append(lines[i][27:33].strip())
                dbus['angV'].append(lines[i][33:40].strip())
                dbus['lMW'].append(lines[i][41:49].strip())
                dbus['lMVAR'].append(lines[i][50:59].strip())
                dbus['gMW'].append(lines[i][60:67].strip())
                dbus['gMVAR'].append(lines[i][68:75].strip())
                dbus['basekV'].append(lines[i][77:83].strip())
                dbus['maxMVAR'].append(lines[i][91:98].strip())
                dbus['minMVAR'].append(lines[i][99:106].strip())
                dbus['shuntG'].append(lines[i][106:114].strip())
                dbus['shuntB'].append(lines[i][115:122].strip())           
                i += 1
                nbus += 1             

        branchTotal = 0      
        if 'BRANCH DATA' in lines[i].strip():    
            i += 1
            while lines[i][0:4].strip() not in endBlock: 
                dbranch['iBus'].append(lines[i][0:4].strip())            
                dbranch['jBus'].append(lines[i][5:9].strip())
                dbranch['area'].append(lines[i][11:12].strip())
                dbranch['type'].append(lines[i][18].strip())
                dbranch['lineR'].append(lines[i][20:29].strip())
                dbranch['lineX'].append(lines[i][30:40].strip())
                dbranch['lineB'].append(lines[i][41:50].strip())            
                dbranch['magTap'].append(lines[i][76:82].strip())            
                dbranch['angTap'].append(lines[i][84:90].strip())        
                i += 1
                branchTotal += 1

    dbus = DF(dbus)
    dbus = dbus.replace(r"^\s*$", '0', regex=True) #transforma vazio em '0'
    dbus = dbus.astype(     
        {
            'num': 'int',
            'name': 'str',
            'area': 'int',
            'type': 'int',
            'magV': 'float',
            'angV': 'float',
            'lMW': 'float',
            'lMVAR': 'float',
            'gMW': 'float',
            'gMVAR': 'float',
            'basekV': 'float',
            'maxMVAR': 'float',
            'minMVAR': 'float', 
            'shuntG': 'float',
            'shuntB': 'float',
        }
    )

    dbranch = DF(dbranch)
    dbranch = dbranch.replace(r"^\s*$", '0', regex=True) #transforma vazio em '0'
    dbranch = dbranch.astype(     
        {
            'iBus': 'int',
            'jBus': 'int',
            'area': 'int',
            'type': 'int',
            'lineR': 'float',
            'lineX': 'float',
            'lineB': 'float',
            'magTap': 'float',
            'angTap': 'float',
        }
    )

    return dbus, dbranch, nbus, MVAbase