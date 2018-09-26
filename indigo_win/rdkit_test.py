from rdkit import Chem
from rdkit.Chem import rdChemReactions
smarts = '[C:1]1[C:2][C:3][C:4]1[C:5]>>[C:1]=[C:2].[C:3]=[C:4][C:5]'
rxn = rdChemReactions.ReactionFromSmarts(smarts)

reacts = (Chem.MolFromSmiles('[C@@H]1(Oc2ccc(C)cc2)[C@H](Oc3ccc(C)cc3)CC1(C)C(=C)C'),)
products = rxn.RunReactants(reacts)

# # S3_1 = "Nc1cc(C(=O)O)c(C(=O)O)cc1S"
# # s3_2 = "Nc1cc(C(=O)O)c(C(=O)O)cc1S.CC=O"
# s1_1 = "N#CC1=CC=C(C2OC2)C=C1"
# s1_3 ="[H]OCC(Cl)c1ccc(C#N)cc1"

for i in range(len(products)):
    print("_____________________________________")
    print(Chem.MolToSmiles(products[i][0]))
    print(Chem.MolToSmiles(products[i][1]))
