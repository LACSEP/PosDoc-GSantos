from numpy import ndarray, zeros
from numpy.linalg import det, eig, solve, inv

def calcYbus(busTotal, dbus, dbranch):  
    #Matriz de admitancias
    yBus: ndarray = zeros(shape=[busTotal,busTotal], dtype='complex')        

    #Linhas e transformadores  
    for _, value in dbranch.iterrows():
        #Elementos fora da diagonal principal       
        iBus = dbus.index[dbus['num'] == value['iBus']][0]
        jBus = dbus.index[dbus['num'] == value['jBus']][0]                       
        if (value['magTap'] == 0.):             
            yBus[iBus,jBus] -= (1/complex(real=value['lineR'],imag=value['lineX'])) 
            yBus[jBus,iBus] -= (1/complex(real=value['lineR'],imag=value['lineX']))                             
            yBus[iBus, iBus] += (1/complex(real=value['lineR'],imag=value['lineX'])) + complex(real=0.,imag=(value['lineB']/2)) 
            yBus[jBus, jBus] += (1/complex(real=value['lineR'],imag=value['lineX'])) + complex(real=0.,imag=(value['lineB']/2))
        else:
            yBus[iBus,jBus] -= (1/complex(real=value['lineR'],imag=value['lineX'])) / float(value['magTap'])
            yBus[jBus,iBus] -= (1/complex(real=value['lineR'],imag=value['lineX'])) / float(value['magTap'])
            yBus[iBus,iBus] += (1/complex(real=value['lineR'],imag=value['lineX'])) / float(value['magTap']) + (1/complex(real=value['lineR'],imag=value['lineX'])) / float(value['magTap']) * (1/float(value['magTap']) -1)
            yBus[jBus,jBus] += (1/complex(real=value['lineR'],imag=value['lineX'])) / float(value['magTap']) + (1/complex(real=value['lineR'],imag=value['lineX'])) * (1 - 1/float(value['magTap']))                
            
    #Banco de capacitores e reatores
    for idx, value in dbus.iterrows():
        if value['shuntB'] != 0.:
            yBus[idx,idx] += complex(0,imag=float(value['shuntB'])) 

    return(yBus)
