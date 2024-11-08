import argparse
from tqdm import tqdm
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize.rdMolStandardize import Normalizer, LargestFragmentChooser, Uncharger

# 对RDKIT的版本有要求

class MolClean(object):
    def __init__(self):
        self.normizer = Normalizer()
        self.lfc = LargestFragmentChooser()
        self.uc = Uncharger()

        # Suppress RDKit warning messages
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
    def clean(self, smi):
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mol = self.normizer.normalize(mol)
            mol = self.lfc.choose(mol)
            mol = self.uc.uncharge(mol)
            smi = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
            return smi
        else:
            return None
 
def main(input_file, output_file):
    mc = MolClean()
 
    with open(input_file, 'r', encoding='utf-8') as f:
        smiles = [r.rstrip() for r in f]
 
    print(f'input SMILES num: {len(smiles)}')
    print('Clean Up Started ...')
 
    mc_smiles = [mc.clean(smi) for smi in tqdm(smiles)]
 
    print('Cleanup Completed')
    print(f'output SMILES num: {len(mc_smiles)}')
 
    with open(output_file, 'w', encoding='utf-8') as f:
        for smi in mc_smiles:
            f.write(smi + '\n')
 
    return
 
if __name__ == '__main__':
    print('cleaner')
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='input file')
    parser.add_argument('output', help='output file')
 
    args = parser.parse_args()
    main(args.input, args.output)
