# smiles-render

A simple CLI that receives smiles and renders its 2D visualization.

## Requirements

- python3

or

- mise

## Setup

```bash
    # If you do not have python3 installed already.
    mise install

    # Config the pip venv.
    python3 -m venv .venv
    pip install -r requirements.txt

    # Config the smiles-render script.
    chmod 751 smiles-render

    # Running the application.
    smiles-render [args]
    
    # or
    python3 src/main [args]
```

## Usage

```bash
    smiles-render -s "NC(O)C(=O)O"
    smiles-render -s "NC(O)C(=O)O" --formmat jpg
    smiles-render -s "NC(O)C(=O)O" -o molecule1.png
    smiles-render -f "molecules.csv" --smiles-column "smiles"
    smiles-render -f "molecules.csv" --smiles-column "smiles" --names-column "molecule name"
    smiles-render -f "molecules.csv" --smiles-column "smiles" --names-column "molecule name" --delimiter ';'
    smiles-render -f "molecules.csv" --smiles-column "smiles" --names-column "molecule name" --delimiter ';' --formmat jpg
```
