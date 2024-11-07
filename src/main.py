from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

m = Chem.MolFromSmiles('NC(O)C(=O)O')

img = Draw.MolToImage(m)
img.save('molecule.png')
