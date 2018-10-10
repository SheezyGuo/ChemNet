import getopt
import sys
from indigo import *


def main(argv):
    reaction = None
    try:
        opts, args = getopt.getopt(argv, "hS:", ["help", "smiles="])
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                print('single_indigo_test.py -S|--smiles')
                sys.exit(0)
            if opt in ("-S", "--smiles"):
                reaction = arg
    except getopt.GetoptError:
        print('Usage: single_indigo_test.py -S|--smiles')
        sys.exit(0)

    indigo = Indigo()
    if not reaction:
        reaction = "c1(cccc(Br)c1)N(=O)=O.c1(ccccc1N)[As](=O)(O)O>>c1(ccccc1[As](=O)(O)O)Nc2cccc(N(=O)=O)c2"
        print("No reaction input,use default reation:{:s}".format(reaction))
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
    main(sys.argv[1:])
