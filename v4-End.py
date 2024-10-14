from time import time
import numpy as np
import multiprocessing as mp

# DMM：E3  coordinate atoms number:6 O-
# FSI: E3  coordinate atoms number:2 O- 3 F- 5 N-

def read_file(file_dir, atoms, box_size):
    with open(file_dir, 'r') as f:
        info = f.readlines() # complet read all data
        # K FSI DMM's :position id mol(special symbol) type(molecules type such as k:type-1)
        # LDMM, LFSI: the RDF of DMM and FSI with K
        posk, idsk = [], []
        molk = {1: {}, 2: {}, 3: {}}
        posDMM, idsDMM, molDMM = [], [], []
        posFSI, idsFSI, molFSI, typeFSI, typeDMM, AGGFSI = [], [], [], {}, {}, {}
        
        # 读取文件中的原子信息
        for i in range(atoms):
            line = info[i].split()
            atom_id, mol, atom_type = int(line[0]), int(line[1]), int(line[2])
            x, y, z = map(float, line[3:6])

            if atom_type == 1:
                posk.append([x, y, z])
                idsk.append(atom_id)
                for key in molk:
                    molk[key].setdefault(str(mol), 0)
            elif atom_type in {6, 8, 9}:
                posDMM.append([x, y, z])
                idsDMM.append(atom_id)
                molDMM.append(mol)
                typeDMM[atom_id] = atom_type
            elif atom_type in {2, 3, 5}:
                posFSI.append([x, y, z])
                idsFSI.append(atom_id)
                molFSI.append(mol)
                typeFSI[atom_id] = atom_type
                AGGFSI.setdefault(str(mol), 0)

    return posk, idsk, molk, posDMM, idsDMM, molDMM, posFSI, idsFSI, molFSI, typeFSI, typeDMM, AGGFSI

def process_file(file, atoms, box_size):
    posk, idsk, molk, posDMM, idsDMM, molDMM, posFSI, idsFSI, molFSI, typeFSI, typeDMM, AGGFSI = read_file(file, atoms, box_size)

    # 转换为数组
    posk = np.array(posk)
    posDMM = np.array(posDMM)
    posFSI = np.array(posFSI)

    total = len(idsk)
    ratioCN, ratioAGGFSI = {}, {}
    ratioAGG = {'AGG': 0, 'CIP': 0, 'SSIP': 0}

    # 遍历K与DMM和FSI的配位数
    for i in range(len(idsk)):
        numDMM, numFSI, numCNDMM, numCNFSI = 0, 0, 0, 0
        molDMM2, molFSI2 = {}, {}

        for j in range(len(idsDMM)):
            xDMM, yDMM, zDMM = posk[i] - posDMM[j]
            xDMM = xDMM - box_size if xDMM > box_size * 0.7 else xDMM + box_size if xDMM <= -box_size * 0.7 else xDMM
            yDMM = yDMM - box_size if yDMM > box_size * 0.7 else yDMM + box_size if yDMM <= -box_size * 0.7 else yDMM
            zDMM = zDMM - box_size if zDMM > box_size * 0.7 else zDMM + box_size if zDMM <= -box_size * 0.7 else zDMM
            LDMM = np.sqrt(xDMM**2 + yDMM**2 + zDMM**2)

            # E3 have 3 Oxygen atoms
            if (typeDMM[idsDMM[j]] == 6 and LDMM <= 4.116):
                if molDMM[j] not in molDMM2:
                    numDMM += 1
                    numCNDMM += 1
                    molDMM2[molDMM[j]] = 'exist'
                else:
                    numCNDMM += 1

        for z in range(len(idsFSI)):
            xFSI, yFSI, zFSI = posk[i] - posFSI[z]
            xFSI = xFSI - box_size if xFSI > box_size * 0.7 else xFSI + box_size if xFSI <= -box_size * 0.7 else xFSI
            yFSI = yFSI - box_size if yFSI > box_size * 0.7 else yFSI + box_size if yFSI <= -box_size * 0.7 else yFSI
            zFSI = zFSI - box_size if zFSI > box_size * 0.7 else zFSI + box_size if zFSI <= -box_size * 0.7 else zFSI

            LFSI = np.sqrt(xFSI**2 + yFSI**2 + zFSI**2)

            # E3 have 3 Oxygen atoms    
            if (typeFSI[idsFSI[z]] == 2 and LFSI <= 4.14) or (typeFSI[idsFSI[z]] == 3 and LFSI <= 3.804) or (typeFSI[idsFSI[z]] == 5 and LFSI <= 6.036):
                if molFSI[z] not in molFSI2:
                    numFSI += 1
                    numCNFSI += 1
                    molFSI2[molFSI[z]] = 'exist'
                else:
                    numCNFSI += 1

        key = f'DMM={numDMM},FSI={numFSI},CNDMM={numCNDMM},CNFSI={numCNFSI}'
        ratioCN[key] = ratioCN.get(key, 0) + 1

    # 遍历FSI与K的配位情况
    for d in range(len(idsFSI)):
        for e in range(len(idsk)):
            xFSI, yFSI, zFSI = posk[e] - posFSI[d]
            xFSI = xFSI - box_size if xFSI > box_size * 0.7 else xFSI + box_size if xFSI <= -box_size * 0.7 else xFSI
            yFSI = yFSI - box_size if yFSI > box_size * 0.7 else yFSI + box_size if yFSI <= -box_size * 0.7 else yFSI
            zFSI = zFSI - box_size if zFSI > box_size * 0.7 else zFSI + box_size if zFSI <= -box_size * 0.7 else zFSI
            LFSI = np.sqrt(xFSI**2 + yFSI**2 + zFSI**2)

            if (typeFSI[idsFSI[d]] == 2 and LFSI <= 4.14) or (typeFSI[idsFSI[d]] == 3 and LFSI <= 3.804) or (typeFSI[idsFSI[d]] == 5 and LFSI <= 6.036): 
                if molk[1].get(str(idsk[e]), 0) == 0:
                    molk[1][str(idsk[e])] = str(molFSI[d])
                    AGGFSI[str(molFSI[d])] += 1
                elif molk[1].get(str(idsk[e])) != str(molFSI[d]):
                    if molk[2].get(str(idsk[e]), 0) == 0:
                        molk[2][str(idsk[e])] = str(molFSI[d])
                        AGGFSI[str(molFSI[d])] += 1
                    elif molk[2].get(str(idsk[e])) != str(molFSI[d]):
                        molk[3][str(idsk[e])] = str(molFSI[d])
                        AGGFSI[str(molFSI[d])] += 1

    # 统计SSIP, CIP, AGG
    for key in molk[1]:
        if molk[1][key] == 0:
            ratioAGG['SSIP'] += 1
        elif molk[2][key] == 0:
            ratioAGG['CIP'] += 1
        else:
            ratioAGG['AGG'] += 1

    for key in AGGFSI:
        agg_key = f'{AGGFSI[key]}K-FSI'
        ratioAGGFSI[agg_key] = ratioAGGFSI.get(agg_key, 0) + 1

    return total, ratioCN, ratioAGG, ratioAGGFSI

# 多进程合并
def collect_results(result, total, ratioCN, ratioAGG, ratioAGGFSI):
    # result: 单个进程结果 
    total_res, ratioCN_res, ratioAGG_res, ratioAGGFSI_res = result
    total.value += total_res
    for key, value in ratioCN_res.items():
        ratioCN[key] = ratioCN.get(key, 0) + value
    for key, value in ratioAGG_res.items():
        ratioAGG[key] = ratioAGG.get(key, 0) + value
    for key, value in ratioAGGFSI_res.items():
        ratioAGGFSI[key] = ratioAGGFSI.get(key, 0) + value

if __name__ == '__main__':
    start = time()
    # 提前定义好参数 file.ipynb会给出
    atoms = 12068
    box_size = 49.99699348576955
    manager = mp.Manager() # 多进程管理器，用于在多进程间共享数据
    # 共享变量 total，赋值0 ，共享ratioCN、ratioAGGFSI、ratioAGG 字典集
    total = manager.Value('i', 0) 
    ratioCN = manager.dict()
    ratioAGGFSI = manager.dict()
    ratioAGG = manager.dict({'AGG': 0, 'CIP': 0, 'SSIP': 0})

    with open('file_list.txt', 'r') as f:
        csv_files = [file.strip() for file in f.readlines()]

    # process pool
    pool = mp.Pool(mp.cpu_count())
    for file in csv_files:
        pool.apply_async(process_file, args=(file, atoms, box_size), callback=lambda res: collect_results(res, total, ratioCN, ratioAGG, ratioAGGFSI))

    pool.close() # close the pool
    pool.join()  # wait all process done

    print('-------end reading files-------')
    print(dict(ratioAGG)) 
    print(dict(ratioAGGFSI)) 
    print(dict(ratioCN)) 

    print('time:', time() - start)
    
    # 将ratioCN 写入文件 如'DMM=3,FSI=0,CNDMM=7,CNFSI=0' 仅写入3 0 7 0
    with open('ratioCN.txt', 'w') as f:
        for key, value in ratioCN.items():
            f.write(f'{key.split(",")[0].split("=")[1]} {key.split(",")[1].split("=")[1]} {key.split(",")[2].split("=")[1]} {key.split(",")[3].split("=")[1]} {value}\n')


