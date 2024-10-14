# Yuan_Scripts
## 如果有好脚本请立即收集并上传到这里， 保持更新 ---- YuanLecun

### 1. preElectro.sh
脚本文件
chkTofchk&predEle.sh 的简化版，只有第(3-5)步骤

### 2. chkTofchk&predEle.sh
脚本文件
(1)遍历本文件夹下所有Guassian计算的chk文件 Loop 
(2)将chk转换为fchk文件
(3)调用Multiwfn_noGUI 对其进行分析，input 12 0 -1 -1 q EOF 生成 out.txt
(4)mini = grep "Minimal value" out.txt | awk -F: '{print $2}' | awk '{print $1}'
(5)maxi = grep "Maximal value" out.txt | awk -F: '{print $3}' | awk '{print $1}'
(6)mini maxi 保存到 result.txt中

### 3. chkTofchk&predEle2.sh
脚本文件：对chkTofchk&predEle.sh存在的问题进行了改进
(1)对chk文件进行了sort，避免了chk文件乱序输出
(2)对所计算文件进行了标注
(3)引入报错机制；在formchk转换失败时，输出错误信息
(4)在Multiwfn_noGUI分析失败时，在result.txt中写入错误文件

### 4. get_out.ipynb
jupyter notebook文件：将out文件夹中最后一次优化后的几何坐标提取出来
(1)遍历当前文件夹下所有的out文件
(1)读取out文件，找到'Standard orientation' 并记录下标
(1)针对每种元素对应的数值 如（C:6 O:8 H:1），进行类别转换并写入为PDB文件

### 5. gjf_get.ipynb
jupyter notebook文件：将PDB文件中的几何坐标提取出来并写入gjf文件
(1)遍历当前文件夹中的PDB文件
(2)逐一读取PDB文件中原子坐标
(3)创建子文件夹存放即将写入的gjf文件
(4)按标准格式写入gjf文件 （几何优化 + 单点能计算）
(5)针对每个gjf文件，按顺序写一个shell脚本文件，用于提交计算

### 6. preddens.sh
脚本文件：利用Multiwfn基于分子表面静电势描述符预测中性分子的晶体密度

### 7. SMI_PDB.ipynb
jupyter notebook文件：读取SMILES文件，利用RDKIT将SMILES转换为分子结构并写为PDB

### 8. Screening.ipynb
jupyter notebook文件：SMI_PDB.ipynb的进阶版
(1)进行了筛选：电荷不为零的分子 和 重原子数大于20的分子 不考虑
(2)仅考虑含N Si P S Cl的分子

### 9. 自动切割表面.ipynb
jupyter notebook文件：接入Materials_Project数据库，并进行晶格表面自动切割
(1)切表面之前要用SpacegroupAnalyzer(struct).get_conventional_standard_structure() 确定已经将晶胞转化成了惯用晶胞
(2)SlabGenerator产生所有可能截断的表面，先读入一个structure类的结构，确定晶面指数：miller_index=[1, 1, 1]，确定slab层的最小厚度：min_slab_size=8.0 (unit: Angstrom)，确定真空层的厚度：min_vacuum_size=15.0 (unit: Angstrom)
(3)在循环所有可能的截断的表面，如果是Au(111)的话只有一个可能的暴露表面，如果是化合物，可能有多种表面
(4)用make_supercell进行扩胞

### 10. get_holumo.sh
脚本文件:调用Multiwfn_noGUI 实现对fchk文件分析, 并提取HOMO和LUMO值存入holumo.txt中
