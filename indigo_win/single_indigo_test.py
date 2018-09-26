from indigo import *


def main():
    indigo = Indigo()
    # rxn = indigo.loadMolecule("ONc1cccc1");
    reaction = "c1(cccc(Br)c1)N(=O)=O.c1(ccccc1N)[As](=O)(O)O>>c1(ccccc1[As](=O)(O)O)Nc2cccc(N(=O)=O)c2"
    print(reaction)
    rxn = indigo.loadReaction(reaction)
    # print("reacting centers:")
    # for m in rxn.iterateMolecules():
    #     for b in m.iterateBonds():
    #         print(rxn.reactingCenter(b))
    # for m in rxn.iterateMolecules():
    #     for b in m.iterateBonds():
    #         rxn.setReactingCenter(b, Indigo.RC_CENTER | Indigo.RC_UNCHANGED)

    loaded_smiles = rxn.smiles()
    print(loaded_smiles)
    rxn.automap("discard")
    automap_smiles = rxn.smiles()
    print(automap_smiles)


if __name__ == "__main__":
    main()
