from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem


def render_smile(smile: str, output_path: str):
    print(f"Rendering {smile} to {output_path}")
    m = Chem.MolFromSmiles(smile)
    img = Draw.MolToImage(m)
    img.save(output_path)
    print("Done!")


def main():
    render_smile("NC(O)C(=O)O", "./molecule.png")


if __name__ == "__main__":
    main()
