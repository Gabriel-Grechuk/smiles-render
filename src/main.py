from rdkit import Chem
from rdkit.Chem import Draw
import argparse
import re


parser = argparse.ArgumentParser(
    prog="smiles-render",
    description="A simple CLI that receives smiles and renders its 2D visualization.",
)
parser.add_argument("-s", "--string", nargs="*", type=str)
parser.add_argument("-o", "--out", nargs="*", type=str)
parser.add_argument("smiles", nargs="*", type=str)
args = parser.parse_args()


def render_smile(smile: str, output_path: str) -> None:
    print(f"Rendering {smile} to {output_path}")
    m = Chem.MolFromSmiles(smile)
    img = Draw.MolToImage(m)
    img.save(output_path)
    print("Done!")


def sanitize_file_name(name: str) -> str:
    return re.sub(r"[^a-zA-Z0-9]", "", name).lower()


def process_args(args: argparse.Namespace) -> None:
    # Positional Argument
    if args.smiles:
        for smile in args.smiles:
            render_smile(smile, f"{sanitize_file_name(smile)}.png")
        return

    # Named output path smiles.
    if len(args.string):
        # Args mismatch.
        if args.out and not len(args.out) == len(args.string):
            raise Exception("Not all smiles input have its correspondent output path")

        # No path smiles.
        elif not args.out:
            for smile in args.string:
                render_smile(smile, f"{sanitize_file_name(smile)}.png")
            return

        # Smiles with path
        else:
            for smile, path in zip(args.string, args.out):
                render_smile(smile, path)
            return


def main():
    process_args(args)


if __name__ == "__main__":
    main()
