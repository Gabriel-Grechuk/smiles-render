from os import read
from rdkit import Chem
from rdkit.Chem import Draw
import argparse
import re
import csv


supported_formats = [
    "BLP",
    "BMP",
    "DDS",
    "DIB",
    "EPS",
    "GIF",
    "ICNS",
    "ICO",
    "IM",
    "JPG",
    "JPEG",
    "MSP",
    "PCX",
    "PFM",
    "PNG",
    "PPM",
    "TIFF",
    "WEBP",
    "XBM",
]


parser = argparse.ArgumentParser(
    prog="smiles-render",
    description="A simple CLI that receives smiles and renders its 2D visualization.",
    usage="""
         %(prog)s -s "NC(O)C(=O)O"
         %(prog)s -s "NC(O)C(=O)O" --format jpg
         %(prog)s -s "NC(O)C(=O)O" -o molecule1.png
         %(prog)s -f "molecules.csv" --smiles-column "smiles"
         %(prog)s -f "molecules.csv" --smiles-column "smiles" --names-column "molecule name"
         %(prog)s -f "molecules.csv" --smiles-column "smiles" --names-column "molecule name" --delimiter ';'
         %(prog)s -f "molecules.csv" --smiles-column "smiles" --names-column "molecule name" --delimiter ';' --format jpg
    """,
)

parser.add_argument(
    "-s", "--string", nargs="*", type=str, help="Smiles that should be rendered."
)

parser.add_argument(
    "-o", "--out", nargs="*", type=str, help="Output path of the rendered smiles."
)

parser.add_argument(
    "-f",
    "--file",
    nargs="*",
    type=str,
    help="The csv file that have the smiles (dependes on --smiles-column --names-column --delimite [optional]",
)

parser.add_argument(
    "--smiles-column",
    nargs="*",
    type=str,
    help="The smiles column in the csv file.",
)

parser.add_argument(
    "--names-column",
    nargs="*",
    type=str,
    help="The molecule or file name column in the csv file for the output file.",
)

parser.add_argument(
    "--delimiter",
    nargs="?",
    type=str,
    default=",",
    help="The csv delimiter token, the default is ','",
)

parser.add_argument(
    "--format",
    nargs="?",
    type=str,
    default="png",
    help="Output format (is overwriten by -o)",
)

parser.add_argument(
    "smiles", nargs="*", type=str, help="The colecules that should be rendered"
)

args = parser.parse_args()


def render_smile(smile: str, output_path: str) -> None:
    print(f"Rendering \033[32m{smile}\033[0m to \033[32m{output_path}\033[0m")
    m = Chem.MolFromSmiles(smile)
    img = Draw.MolToImage(m)
    img.save(output_path)
    print("Done!")


def sanitize_file_name(name: str) -> str:
    return re.sub(r"[^a-zA-Z0-9]", "", name).lower()


def read_file(
    path: str, delimiter: str, smiles_column: str, names_column: str | None = None
) -> list[tuple[str, str]]:
    file = open(path, newline="")
    csv_data = csv.DictReader(file, delimiter=delimiter)

    if not csv_data:
        raise Exception("Could not read csv file")

    if csv_data.fieldnames is not None:
        if smiles_column not in csv_data.fieldnames:
            raise Exception(f"{smiles_column} column not found")
        if names_column is not None:
            if names_column not in csv_data.fieldnames:
                raise Exception(f"{names_column} column not found")
    else:
        raise Exception("Could not read csv file")

    smiles: list[str] = []
    names: list[str] = []

    for row in csv_data:
        smiles.append(row[smiles_column].strip())
        if names_column:
            names.append(row[names_column].strip())
        else:
            names.append("")

    return list(zip(smiles, names))


def process_args(args: argparse.Namespace) -> None:
    # Check the selected output format.
    if args.format:
        if args.format.upper() not in supported_formats:
            raise Exception(f"Image format {args.format.upper()} not supported")

    # Positional Argument
    if args.smiles:
        for smile in args.smiles:
            render_smile(smile, f"{sanitize_file_name(smile)}.{args.format}")

        return

    # Named output path smiles.
    if args.string:
        # Args mismatch.
        if args.out and not len(args.out) == len(args.string):
            raise Exception("Not all smiles input have its correspondent output path")

        # No path smiles.
        elif not args.out:
            for smile in args.string:
                render_smile(smile, f"{sanitize_file_name(smile)}.{args.format}")

            return

        # Smiles with path.
        else:
            for smile, path in zip(args.string, args.out):
                render_smile(smile, path)

            return

    # Generate from cvs file.
    if args.file:
        if not args.smiles_column:
            raise Exception("Argument --smiles-column missing")

        data = read_file(
            path=args.file[0],
            smiles_column=args.smiles_column[0],
            delimiter=args.delimiter[0],
            names_column=args.names_column[0] if args.names_column else None,
        )

        for smile, name in data:
            if name:
                render_smile(smile, f"{sanitize_file_name(name)}.{args.format}")
            else:
                render_smile(smile, f"{sanitize_file_name(smile)}.{args.format}")

        return

    parser.print_help()


def main():
    process_args(args)


if __name__ == "__main__":
    main()
