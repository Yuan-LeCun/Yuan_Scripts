import numpy as np
import os
from time import time

start = time()

criterion = []
line_num, numtotal, atoms = 0, 0, 12068
boxsizex, boxsizey, boxsizez = 49.99699348576955, 49.9806972036557, 49.94810463946648

with open('ratioCN.txt', 'r') as file:
    for line in file:   line_num += 1
    file = 'ratioCN.txt'
    if os.path.exists('structure'):
        pass
    else: os.mkdir('structure')
    
with open(file, 'r') as file:
    info = file.readlines()
    for i in range(line_num):
        line = info[i].split()
        p1, p2, p3, p4 = int(line[0]), int(line[1]), int(line[2]), int(line[3])
        criterion.append([p1, p2, p3, p4])
    criterion = np.array(criterion)

def read_file(file_dir):
    with open(file_dir, 'r') as f:

        posK, IdsK, molK = [], [], []
        posEC, IdsEC, molEC, typeEC = [], [], [], []
        posFSI, IdsFSI, molFSI, typeFSI = [], [], [], []
        
        for _ in range(atoms):
            Id, Mol, Type, x, y, z = map(float, f.readline().split())
            if Type == 1: 
                posK.append([x, y, z])
                IdsK.append(Id)
                molK.append(Mol)
            elif Type in {6, 7, 8, 9, 10}:
                posEC.append([x, y, z])
                IdsEC.append(Id)
                typeEC.append(Type)
                molEC.append(Mol)
            elif Type in {2, 3, 4, 5}:
                posFSI.append([x, y, z])
                IdsFSI.append(Id)
                typeFSI.append(Type)
                molFSI.append(Mol)
        posK, posEC, posFSI = np.array(posK), np.array(posEC), np.array(posFSI)
    return posK, IdsK, molK, posEC, IdsEC, molEC, typeEC, posFSI, IdsFSI, molFSI, typeFSI

def process_file(file_dir, num, count, cycle):
    posK, IdsK, molK, posEC, IdsEC, molEC, typeEC, posFSI, IdsFSI, molFSI, typeFSI = read_file(file_dir)
    for a in range(len(IdsK)):
        numEC, numFSI, numCN = 0, 0, 0
        molDEC2, molEC2, molFSI2 = {}, {}, {}
        MOL_DEC, MOL_EC, MOL_FSI = [], [], []

        for c in range(len(IdsEC)):
            if criterion[num, 0] == 0:    
                break
            xEC, yEC, zEC = posK[a] - posEC[c]
            xEC = xEC - boxsizex if xEC > boxsizex * 0.7 else xEC + boxsizex if xEC < boxsizex * -0.7 else xEC
            yEC = yEC - boxsizey if yEC > boxsizey * 0.7 else yEC + boxsizey if yEC < boxsizey * -0.7 else yEC
            zEC = zEC - boxsizez if zEC > boxsizez * 0.7 else zEC + boxsizez if zEC < boxsizez * -0.7 else zEC
            LEC = np.sqrt(xEC**2 + yEC**2 + zEC**2)

            # EC : E3
            if typeEC[c] == 6 and LEC <= 4.116:
                if molEC[c] not in molEC2:
                    numEC += 1
                    numCN += 1
                    molEC2[molEC[c]] = 'exist'
                    MOL_EC.append(molEC[c])
                else:
                    numCN += 1
        if numEC != criterion[num, 0]:
            continue                      

        for d in range(len(IdsFSI)):
            if criterion[num, 1] == 0:
                break

            xFSI, yFSI, zFSI = posK[a] - posFSI[d]
            xFSI = xFSI - boxsizex if xFSI > boxsizex * 0.7 else xFSI + boxsizex if xFSI < boxsizex * -0.7 else xFSI
            yFSI = yFSI - boxsizey if yFSI > boxsizey * 0.7 else yFSI + boxsizey if yFSI < boxsizey * -0.7 else yFSI
            zFSI = zFSI - boxsizez if zFSI > boxsizez * 0.7 else zFSI + boxsizez if zFSI < boxsizez * -0.7 else zFSI
            LFSI = np.sqrt(xFSI**2 + yFSI**2 + zFSI**2)
            if (typeFSI[d] == 2 and LFSI <= 4.14) or (typeFSI[d] == 3 and LFSI <= 3.804) or (typeFSI[d] == 5 and LFSI <= 6.036 ):
                if molFSI[d] not in molFSI2:
                    numFSI += 1
                    numCN += 1
                    molFSI2[molFSI[d]] = 'exist'
                    MOL_FSI.append(molFSI[d])
                else:
                    numCN += 1
            
        if numFSI != criterion[num,1]:
            continue
        
        else:       
            count += 1
            global numtotal
            numtotal += 1
            os.chdir('./structure')
            f = open(f'str{numtotal}.gjf', 'w')
            f.write(f'%nprocshared=48\n%mem=80GB\n%chk=\\str{count}.chk\n')
            f.write('# B3LYP/TZVP em=GD3BJ scrf=(iefpcm,read) opt \n\n')
            f.write('Title Card Required\n\n')
            charge = 1 - criterion[num,1]
            f.write(f'{charge} 1\n K          {posK[a,0]}   {posK[a,1]}   {posK[a,2]}\n')

            if criterion[num,0] != 0:
                for item in MOL_EC:
                    count_EC = 0
                    for c in range(len(IdsEC)):
                        if count_EC < 10 and item == molEC[c]:
                            count_EC += 1
                            f.write(' O' if typeEC[c] == 6 else ' C' if typeEC[c] in {7, 8} else ' H')
                            x = posEC[c,0] + boxsizex * (1 if posK[a,0] - posEC[c,0] > boxsizex * 0.7 else -1 if posK[a,0] - posEC[c,0] < -boxsizex * 0.7 else 0)
                            y = posEC[c,1] + boxsizey * (1 if posK[a,1] - posEC[c,1] > boxsizey * 0.7 else -1 if posK[a,1] - posEC[c,1] < -boxsizey * 0.7 else 0)
                            z = posEC[c,2] + boxsizez * (1 if posK[a,2] - posEC[c,2] > boxsizez * 0.7 else -1 if posK[a,2] - posEC[c,2] < -boxsizez * 0.7 else 0)
                            f.write(f'          {x}   {y}   {z}\n')
            
            if criterion[num,1] != 0:
                for items in MOL_FSI:
                    count_FSI = 0
                    for d in range(len(IdsFSI)):
                        if count_FSI < 9 and items == molFSI[d]:
                            count_FSI += 1
                            f.write(' O' if typeFSI[d] == 2 else ' F' if typeFSI[d] == 3 else ' S' if typeFSI[d] == 4 else ' N')
                            x = posFSI[d,0] + boxsizex * (1 if posK[a,0] - posFSI[d,0] > boxsizex * 0.7 else -1 if posK[a,0] - posFSI[d,0] < -boxsizex * 0.7 else 0)
                            y = posFSI[d,1] + boxsizey * (1 if posK[a,1] - posFSI[d,1] > boxsizey * 0.7 else -1 if posK[a,1] - posFSI[d,1] < -boxsizey * 0.7 else 0)
                            z = posFSI[d,2] + boxsizez * (1 if posK[a,2] - posFSI[d,2] > boxsizez * 0.7 else -1 if posK[a,2] - posFSI[d,2] < -boxsizez * 0.7 else 0)
                            f.write(f'          {x}   {y}   {z}\n')
            f.write(f'\nsurface=sas\neps=69.43398271\nepsinf=1.886804914\n\n--link1--\n%nprocshared=48\n%mem=80GB\n%oldchk=\\str{count}.chk\n%chk=\\str{count}_2.chk\n# B3LYP/gen em=GD3BJ scrf=(iefpcm,read) geom=allcheck guess=read\n\n@./ma-TZVP.txt/N\n\neps=69.43398271\nepsinf=1.886804914\n\n')
            f.close()
            os.chdir('../')
            print('%d scanned' %cycle)
            cycle += 1
            break   
    return count, cycle



if __name__ == '__main__':

    for num in range(line_num):
        count, countmax, cycle = 0, 1, 1

        with open('file_list.txt', 'r') as f:
            csv_files = f.readlines()
            csv_files = [csv_file.strip() for csv_file in csv_files]
            
        for file in csv_files:
            count, cycle = process_file(file, num, count, cycle)
            if count >= countmax:
                break
        print('%d criterion completed' %(num+1))
    print('all completed')
    print(f'Time taken: {time() - start} seconds')