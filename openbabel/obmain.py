# -*- coding: utf-8 -*-
"""
__title__ = 'openbabel'
__package = ''
__project__ ='ChemNet"
__author__ = 'gsh'
__date_ = '2018.10.26'
__IDE__ = 'PyCharm'
# code is far away from bugs with the god animal protecting
    I love animals. They taste delicious.
             ┏┓   ┏┓
            ┏┛┻━━━┛┻┓
            ┃       ┃
            ┃┳┛  ┗┳ ┃
            ┃┻    ┻ ┃
            ┗━┓   ┏━┛
              ┃   ┗━━━┓
              ┃神兽保佑 ┣┓
              ┃永无BUG ┏┛
              ┗┓┓┏━┳┓┏┛
               ┃┫┫ ┃┫┫
               ┗┻┛ ┗┻┛
"""

import numpy as np
import openbabel as ob


def get_adjacent_matrix(smiles):
    obConversion = ob.OBConversion()
    mol = ob.OBMol()
    obConversion.SetInFormat("smi")
    obConversion.ReadString(mol, smiles)
    M = np.zeros([mol.NumAtoms(), mol.NumAtoms()])
    for bond in ob.OBMolBondIter(mol):
        start = bond.GetBeginAtom().GetId()
        end = bond.GetEndAtom().GetId()
        order = bond.GetBondOrder()
        M[start, end] = order
        M[end, start] = order
    return M


def get_adjacent_matrix_by_sdf(smiles):
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("smi", "sdf")
    mol = ob.OBMol()
    obConversion.ReadString(mol, smiles)
    sdfString = obConversion.WriteString(mol)
    obConversion.SetInFormat("sdf")
    obConversion.ReadString(mol, sdfString)
    M = np.zeros([mol.NumAtoms(), mol.NumAtoms()])
    for bond in ob.OBMolBondIter(mol):
        start = bond.GetBeginAtom().GetId()
        end = bond.GetEndAtom().GetId()
        order = bond.GetBondOrder()
        M[start, end] = order
        M[end, start] = order
    return sdfString, M


def get_neighbor_graph(matrix, index, dist=3):
    assert type(matrix) == np.ndarray
    if index > matrix.shape[0] or index < 0:
        raise ValueError(index)
    neighbors = dict()
    neighbors["neighbor_1"] = set()
    neighbors["neighbor_1"].add(index)
    index_list = np.array(range(M.shape[0]))
    for i in range(2, dist + 1):
        prev_key = "neighbor_{:d}".format(i - 1)
        cur_key = "neighbor_{:d}".format(i)
        neighbors[cur_key] = set()
        prev_points = neighbors[prev_key]
        for point in prev_points:
            for new_point in index_list[M[point, :] > 0].tolist():
                neighbors[cur_key].add(new_point)
    target_points = set()
    for point_set in neighbors.values():
        for point in point_set:
            target_points.add(point)
    target_points = list(target_points)
    return M[target_points, :][:, target_points]


if __name__ == "__main__":
    # smiles = "C1(C)(C)c2c(cccc2)N(c3ccc(N(c4ccc(N5c6c(cccc6)C(C)(C)c(cccc7)c57)cc4)c8ccc([nH]c(ccc(c9ccccc9)c%10)c%10%11)c%11c8)cc3)c(cccc%12)c1%12"
    smiles = "c12c(ccc(Br)c1)C=Cc(cccc3)c3C2=O"
    sdf, M = get_adjacent_matrix_by_sdf(smiles)
    print(smiles)
    print(sdf)
    print(M)
    print(get_neighbor_graph(M, 2))
