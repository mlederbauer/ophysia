import argparse
from ophysia.example_module import hello_smiles

def main() -> None:
    parser = argparse.ArgumentParser(description='Process some smiles.')
    parser.add_argument('smiles', type=str, help='A text string representing a SMILES notation or any string.')

    args = parser.parse_args()

    result = hello_smiles(args.smiles)
    print(result)

if __name__ == "__main__":
    main()
