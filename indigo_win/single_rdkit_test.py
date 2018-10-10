import sys
from rdkit import Chem

import getopt
from rdkit.Chem import rdChemReactions


def main(argv):
    smiles = None
    smarts = None
    try:
        opts, args = getopt.getopt(argv, "L:T:", ["smarts=", "smiles="])
        for opt, arg in opts:
            if opt in ["-L", "--smiles"]:
                smiles = arg
            if opt in ["-T", "--smarts"]:
                smarts = arg
    except getopt.GetoptError:
        print("Usage:python single_rdkit_test.py -L|--smiles smiles -T|--smarts smarts")
    if not smiles:
        smiles = '[C@@H]1(Oc2ccc(C)cc2)[C@H](Oc3ccc(C)cc3)CC1(C)C(=C)C'
        print("Use default smiles:{:s}".format(smiles))
    if not smarts:
        smarts = '[C:1]1[C:2][C:3][C:4]1[C:5]>>[C:1]=[C:2].[C:3]=[C:4][C:5]'
        print("Use default smarts:{:s}".format(smarts))
    rxn = rdChemReactions.ReactionFromSmarts(smarts)
    reacts = (Chem.MolFromSmiles(smiles),)
    products = rxn.RunReactants(reacts)

    # # S3_1 = "Nc1cc(C(=O)O)c(C(=O)O)cc1S"
    # # s3_2 = "Nc1cc(C(=O)O)c(C(=O)O)cc1S.CC=O"
    # s1_1 = "N#CC1=CC=C(C2OC2)C=C1"
    # s1_3 ="[H]OCC(Cl)c1ccc(C#N)cc1"

    for i in range(len(products)):
        print("_____________________________________")
        print(Chem.MolToSmiles(products[i][0]))
        print(Chem.MolToSmiles(products[i][1]))


if __name__ == "__main__":
    main(sys.argv[1:])
