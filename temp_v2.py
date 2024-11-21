import numpy as np
import time

"""
溶剂化壳层中离子键平均寿命!!!
该代码用于跟踪钾离子与周围分子(H2O、DMSO、FSI)之间的连接关系
通过比较当前和之前的连接分子集合，检测连接的变化，并计算连接的持续时间，最终求得平均连接寿命
"""

def calculate_distance(pos1, pos2, boxsize):
    delta = pos1 - pos2
    delta = np.where(delta > boxsize * 0.7, delta - boxsize, delta)
    delta = np.where(delta < -boxsize * 0.7, delta + boxsize, delta)
    return np.linalg.norm(delta, axis=1)

def read_file(filename, atoms):
    with open(filename, 'r') as f:
        info = f.readlines()
        posK, IdsK = [], []
        posH2O, IdsH2O = [], []
        posDMSO, IdsDMSO = [], []
        posFSI, IdsFSI, molFSI, typeFSI = [], [], [], []

        for i in range(atoms):
            line = info[i].split()
            Id, Mol, Type = int(line[0]), int(line[1]), int(line[2])
            x, y, z = map(float, line[3:6])
            if Type == 7:
                posK.append([x, y, z])
                IdsK.append(Id)
            elif Type == 1:
                posH2O.append([x, y, z])
                IdsH2O.append(Id)
            elif Type == 5:
                posDMSO.append([x, y, z])
                IdsDMSO.append(Id)
            elif Type in {8, 9, 11}:
                posFSI.append([x, y, z])
                IdsFSI.append(Id)
                typeFSI.append(Type)
                molFSI.append(Mol)

    return np.array(posK), IdsK, np.array(posH2O), IdsH2O, np.array(posDMSO), IdsDMSO, np.array(posFSI), IdsFSI, molFSI, np.array(typeFSI)

start = time.time()

def main():
    total, cycle = 0, 1
    K_all, K_H2O, K_DMSO, K_FSI= {}, {}, {}, {}
    time_all, time_H2O, time_DMSO, time_FSI = 0, 0, 0, 0
    N_all, N_H2O, N_DMSO, N_FSI = 0, 0, 0, 0
    boxsize = 52.68883463461055

    with open('file_list.txt', 'r') as f:
        csv_files = f.readlines()
    for file in csv_files:
        file_dir = file.strip()
        posK, IdsK, posH2O, IdsH2O, posDMSO, IdsDMSO, posFSI, IdsFSI, molFSI, typeFSI = read_file(file_dir, 11965)
        total += len(IdsK)
        
        step = int(file_dir.split('.')[-2])
        time_step = int(file_dir.split('.')[-2])/500

        # 特殊判断: 步数为0时肯定键的平均寿命也为0
        if step == 0:
            # IdsK是所有钾离子的ID列表
            # 对于每个钾离子，创建一个包含当前时间步长 time_step 的列表，用于记录连接信息
            for IdK in IdsK:
                K_all[IdK] = [time_step]
                K_H2O[IdK] = [time_step]
                K_DMSO[IdK] = [time_step]
                K_FSI[IdK] = [time_step]
        # enumerate() get index and value
        for idx, IdK in enumerate(IdsK):
            pos_K = posK[idx]
            # 计算距离
            dist_H2O = calculate_distance(posH2O, pos_K, boxsize)
            dist_DMSO = calculate_distance(posDMSO, pos_K, boxsize)
            dist_FSI = calculate_distance(posFSI, pos_K, boxsize)

            # 距离筛选 根据设定的阈值和相吸引的atom_type来获取ID列表
            connected_H2O = set(np.array(IdsH2O)[dist_H2O <= 3.564])
            connected_DMSO = set(np.array(IdsDMSO)[dist_DMSO <= 3.684])
            cutoff_FSI = np.where(
                typeFSI == 8, 3.684,
                np.where(
                    typeFSI == 9, 3.516,
                    3.66    # 除了8,9以外的数值都用3.66作为阈值
                )
            )
            connected_FSI = set(np.array(IdsFSI)[dist_FSI <= cutoff_FSI])
            # 合体!
            connected_all = connected_H2O | connected_DMSO | connected_FSI
            # print('prev_all:', prev_all)
            # print('connected_all:', connected_all)
            # 更新生命期
            if step != 0:
                # 获取之前连接的分子的集合(判断是否存在)
                prev_all = set(K_all[IdK][1:]) if len(K_all[IdK]) > 1 else set()
                # print('prev_all:', prev_all)
                prev_H2O = set(K_H2O[IdK][1:]) if len(K_H2O[IdK]) > 1 else set()
                prev_DMSO = set(K_DMSO[IdK][1:]) if len(K_DMSO[IdK]) > 1 else set()
                prev_FSI = set(K_FSI[IdK][1:]) if len(K_FSI[IdK]) > 1 else set()
                # 比较当前连接的分子集合 与之前有什么不同?
                # 若不同, 表明连接发生了变化, 更新生命期和连接信息
                # K_all[Idk][0]表示当前分子的开始连接时间, 那么time_step - K_all[Idk][0]就是当前分子的连接持续时间
                if prev_all.isdisjoint(connected_all):
                    N_all += 1          # 表示持续了几次?
                    time_all += (time_step - K_all[IdK][0])
                    K_all[IdK] = [time_step] + list(connected_all)
                if prev_H2O != connected_H2O:
                    N_H2O += 1
                    time_H2O += (time_step - K_H2O[IdK][0])
                    K_H2O[IdK] = [time_step] + list(connected_H2O)
                if prev_DMSO != connected_DMSO:
                    N_DMSO += 1
                    time_DMSO += (time_step - K_DMSO[IdK][0])
                    K_DMSO[IdK] = [time_step] + list(connected_DMSO)
                if prev_FSI != connected_FSI:
                    N_FSI += 1
                    time_FSI += (time_step - K_FSI[IdK][0])
                    K_FSI[IdK] = [time_step] + list(connected_FSI)
            else:
                K_all[IdK] += list(connected_all)
                K_H2O[IdK] += list(connected_H2O)
                K_DMSO[IdK] += list(connected_DMSO)
                K_FSI[IdK] += list(connected_FSI)
        print(f'{cycle} completed')
        # print(K_all[IdK][0])
        cycle += 1

    stop = time.time()
    # 计算平均寿命
    ave_all = time_all / N_all
    ave_H2O = time_H2O / N_H2O if N_H2O else 0
    ave_DMSO = time_DMSO / N_DMSO if N_DMSO else 0
    ave_FSI = time_FSI / N_FSI if N_FSI else 0

    print(f'life time of all is {ave_all} ps')
    print(f'life time of H2O is {ave_H2O} ps')
    print(f'life time of DMSO is {ave_DMSO} ps')
    print(f'life time of FSI is {ave_FSI} ps')
    print(f'total time: {stop - start} s')
    print('time all:', time_all)
    print('N all:', N_all)
    # print(N_H2O + N_DMSO + N_FSI == N_all)
    # print(N_DMSO, N_H2O, N_FSI)
    print('time DMSO:', time_DMSO)
    print('N_DMSO:', N_DMSO)
    print('N_H2O:', N_H2O)
    print('N_FSI:', N_FSI)
if __name__ == '__main__':
    main()
