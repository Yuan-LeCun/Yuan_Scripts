# Yuan_Scripts
## 如果有好脚本请立即收集并上传到这里， 保持更新 ---- YuanLecun

### 1. SMI_PDB.ipynb
jupyter notebook文件：读取SMILES文件，利用RDKIT将SMILES转换为分子结构并写为PDB

### 2. Screening.ipynb
jupyter notebook文件：SMI_PDB.ipynb的进阶版
(1)进行了筛选：电荷不为零的分子 和 重原子数大于20的分子 不考虑
(2)仅考虑含N Si P S Cl的分子

### 3. cas_smiles.ipynb
调用Pubchem API直接将cas转换成smiles

### 4. PDB_SMI.ipynb
jupyter notebook文件:读取所需PDB文件(Screening.txt目录),创建result.txt 逐行读取PDB文件内容，将其转换为SMILES描述符并写入result.txt中

### 5. gjf_get.ipynb
jupyter notebook文件：将PDB文件中的几何坐标提取出来并写入gjf文件
(1)遍历当前文件夹中的PDB文件
(2)逐一读取PDB文件中原子坐标
(3)创建子文件夹存放即将写入的gjf文件
(4)按标准格式写入gjf文件 （几何优化 + 单点能计算）
(5)针对每个gjf文件，按顺序写一个shell脚本文件，用于提交计算

### 6. smi_gjf_get.ipynb
### 7. smi_gjf_get_local.ipynb
本质上是一样的, 只是一个是为本地写gjf文件, 一个是为服务器上写gjf文件, 读取smiles文件时形式偶尔需要修改!

### 8. get_out.ipynb
jupyter notebook文件：将out文件夹中最后一次优化后的几何坐标提取出来
(1)遍历当前文件夹下所有的out文件
(1)读取out文件，找到'Standard orientation' 并记录下标
(1)针对每种元素对应的数值 如（C:6 O:8 H:1），进行类别转换并写入为PDB文件

### 9. get_wrongoutPDB_gjf.py
get_out.py的进阶版
(1)读取所有out文件, 报错的out文件则拿出其最后一次几何优化后的原子坐标信息并生成PDB
(2)需要注意的是在写成PDB时,需要注意原子类型的添加
(3)最后写成gjf, bash, gaussslurm.sh

### 10. preElectro.sh
脚本文件
chkTofchk&predEle.sh 的简化版，只有第(3-5)步骤

### 11. chkTofchk&predEle.sh
脚本文件
(1)遍历本文件夹下所有Guassian计算的chk文件 Loop 
(2)将chk转换为fchk文件
(3)调用Multiwfn_noGUI 对其进行分析，input 12 0 -1 -1 q EOF 生成 out.txt
(4)mini = grep "Minimal value" out.txt | awk -F: '{print $2}' | awk '{print $1}'
(5)maxi = grep "Maximal value" out.txt | awk -F: '{print $3}' | awk '{print $1}'
(6)mini maxi 保存到 result.txt中

#### 11.1. chkTofchk&predEle2.sh
脚本文件：对chkTofchk&predEle.sh存在的问题进行了改进
(1)对chk文件进行了sort，避免了chk文件乱序输出
(2)对所计算文件进行了标注
(3)引入报错机制；在formchk转换失败时，输出错误信息
(4)在Multiwfn_noGUI分析失败时，在result.txt中写入错误文件

#### 11.2 chkTofchk&predEle3.sh
对chkTofchk&predEle2.sh进一步修改, 加入了MPI, HOMO, LUMO

#### 11.3 chkTofchk&predEle5.sh
对chkTofchk&predEle4.sh进一步修改, 加入了energy

### 12. preddens.sh
脚本文件：利用Multiwfn基于分子表面静电势描述符预测中性分子的晶体密度

### 13. get_holumo.sh
脚本文件:调用Multiwfn_noGUI 实现对本文件夹下的所有fchk文件分析, 并提取HOMO和LUMO值存入holumo.txt中

### 14. 自动切割表面.ipynb
jupyter notebook文件：接入Materials_Project数据库，并进行晶格表面自动切割
(1)切表面之前要用SpacegroupAnalyzer(struct).get_conventional_standard_structure() 确定已经将晶胞转化成了惯用晶胞
(2)SlabGenerator产生所有可能截断的表面，先读入一个structure类的结构，确定晶面指数：miller_index=[1, 1, 1]，确定slab层的最小厚度：min_slab_size=8.0 (unit: Angstrom)，确定真空层的厚度：min_vacuum_size=15.0 (unit: Angstrom)
(3)在循环所有可能的截断的表面，如果是Au(111)的话只有一个可能的暴露表面，如果是化合物，可能有多种表面
(4)用make_supercell进行扩胞

### 15. file.ipynb
数据提取,并写入atom.{step}.csv文件和file_list.txt

### 15.1.file_v2.ipynb
file的进阶版, 写成函数方便后续操作

### 16.v4-End.py 
全CPU加速,实现溶剂化壳层分类与统计
(1)在MD后实现了阴阳离子和溶剂分子的划分,如放入posk posDMM posFSI中
(2)两个for嵌套,实现碱金属离子与溶剂分子 和 阴离子的配位统计(溶剂化壳层内,RDF距离限制)
(3)两个for嵌套,实现阴离子与碱金属离子的配位统计
(4)统计ratioAGG:SSIP,CIP,AGG(FSI-*K)
(5)统计ratioAGGFSI:*FSI-K
(6)统计ratioCN(such as DMM=3 FSI=0 CNDMM=7 CNFSI=0) 溶剂化壳层内溶剂分子数量及溶剂分子上配位数

#### 16.1 v5.py
v4-end.py进阶版
对变量名统一进行了规范, 如cation, anion, solvent_e3, solvent_n7, solvent_h6

#### 16.2 v5_single_solvent.py
变量名规范后, 主要针对单溶剂

#### 16.3 v5_three_solvent.py
针对三溶剂

### 17.str_outputs_3.py
根据ratioCN:配位情况及配位数 将对应的溶剂化团簇提取出来生成gjf文件

### 18. smiles&graph.ipynb
(1)利用NetworkX 和 RDKit包读取SMILES
(2)结构绘制为图片并保存
(3)SMILES与图片一同插入进excel中
(4)第三列, 第四列分别是分子式和分子重量用于描述分子

### 19. answer&que.ipynb
(1)与smiles&graph.ipynb联用, 读取dataset_with_image.xlsx文件
(2)通过提示栏输入分子式和分子重量, 将对应的Smiles文件提取进extracted_smiles.xlsx文件中, 用于后续的Smiles_trans_PDB_trans_gjf

### 20. Cleaner.py
(1)对RDKIT版本有特殊要求, 可能会报错!
(2)将ISOsmiles转换成CanonicalSmiles, 对Smiles进行清洗和标准化格式转换

### 21.File_CarbonOxgen_addLi.ipynb
(1) 自动识别PDB并提取PDB原子坐标信息和原子类别
(2)在羰基氧附近添加一个Li原子, 且保持Li原子与羰基氧的距离为3埃, Li原子远离其他原子

### 22.MultiwfnMLhelper.py
Notice : 该py文件仅在Linux环境下使用, 且需要能直接调用Multiwfn or Multiwfn_noGUI , setting.ini中i_silence设置为1
(1)调用Multiwfn对fchk文件进行解析, 并提取相应的分子性质描述符
(2) 众多分子性质描述符可用于后续的机器学习

### 23.anion_v4.py
Py文件, 计算以FSI为中心的溶剂化壳层周围的成份(FSI-K, FSI-FSI, FSI-Solvent)

#### 23.1 anion_v6.py
anion_v4.py进阶版, 全CPU 加速!

### 24. temp_v2.py
py文件, 计算溶剂化壳层中离子键平均寿命!
该代码用于跟踪钾离子与周围分子(H2O、DMSO、FSI)之间的连接关系
通过比较当前和之前的连接分子集合，检测连接的变化，并计算连接的持续时间，最终求得平均连接寿命

#### 24.1 temp_v3.py
temp_v2.py进阶版, 全CPU加速




