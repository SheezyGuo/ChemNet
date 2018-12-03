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

import re
import sys

import numpy as np
import openbabel as ob
import pybel


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
    return M, mol


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
    return M, sdfString


def get_neighbor_graph(matrix, index, dist=3):
    assert type(matrix) == np.ndarray
    if index > matrix.shape[0] or index < 0:
        raise ValueError(index)
    neighbors = dict()
    neighbors["neighbor_1"] = set()
    neighbors["neighbor_1"].add(index)
    index_list = np.array(range(matrix.shape[0]))
    for i in range(2, dist + 1):
        prev_key = "neighbor_{:d}".format(i - 1)
        cur_key = "neighbor_{:d}".format(i)
        neighbors[cur_key] = set()
        prev_points = neighbors[prev_key]
        for point in prev_points:
            for new_point in index_list[matrix[point, :] > 0].tolist():
                neighbors[cur_key].add(new_point)
    target_points = set()
    for point_set in neighbors.values():
        for point in point_set:
            target_points.add(point)
    target_points = list(target_points)
    target_points.sort()
    return matrix[target_points, :][:, target_points], target_points


def get_sub_smiles(smiles, points=None, pos=None):
    M, mol = get_adjacent_matrix(smiles)
    print("分子内部顺序")
    for atom in ob.OBMolAtomIter(mol):
        print(atom.GetId(), atom.GetIdx(), atom.GetIndex(), atom.GetAtomicNum(), atom.GetType())
    if points is None and pos is not None:
        sub_M, points = get_neighbor_graph(M, pos)
    if isinstance(points, set):
        points = list(points)
    points.sort()
    new_mol = ob.OBMol()
    for p in points:
        atom = mol.GetAtomById(p)
        # atom.
        new_mol.AddAtom(atom)
    for i in range(len(points)):
        for j in range(len(points)):
            if M[points[i], points[j]] > 0:
                new_mol.AddBond(i + 1, j + 1, int(M[points[i], points[j]]))
    OBConversion = ob.OBConversion()
    print("中心分子内部顺序")
    for atom in ob.OBMolAtomIter(new_mol):
        print(atom.GetId(), atom.GetIdx(), atom.GetIndex(), atom.GetAtomicNum(), atom.GetType())
    OBConversion.SetOutFormat("smi")
    new_mol_smiles = OBConversion.WriteString(new_mol).strip()
    print("中心smiles")
    print(new_mol_smiles)
    print("#" * 80)
    return new_mol_smiles


def get_index_map_from_smarts(mapped_smarts):
    assert len(mapped_smarts.split(">>")) == 2
    # left hand side   right hand side
    parts = mapped_smarts.split(">>")
    lhs = parts[0]
    rhs = parts[1]
    reactants = lhs.split(".")
    products = rhs.split(".")
    pattern_str = r"\[[A-Z0-9]+:(\d+)\]|(Mg|Br|[A-Z]+[a-z]*)"
    pattern = re.compile(pattern_str)
    get_record_index = lambda i_, n_: "{:d}_{:d}".format(i_, n_)
    reactant_map = dict()
    product_map = dict()
    for i, reactant in enumerate(reactants):
        count = 0
        iter = re.finditer(pattern, reactant)
        for matcher in iter:
            automap_index = matcher.group(1)
            if automap_index:
                reactant_map[str(automap_index)] = get_record_index(i, count)
            else:
                print("Atom not matched:{:s},pos:{:d},in {:s}".format(matcher.group(2), count, reactant),
                      file=sys.stderr)
            count = count + 1
    for i, product in enumerate(products):
        count = 0
        iter = re.finditer(pattern, product)
        for matcher in iter:
            automap_index = matcher.group(1)
            if automap_index:
                product_map[str(automap_index)] = get_record_index(i, count)
            else:
                print("Atom not matched:{:s},pos:{:d},in {:s}".format(matcher.group(2), count, product),
                      file=sys.stderr)
            count = count + 1
    reverse_reactant_map = dict()
    reverse_product_map = dict()
    for key, value in reactant_map.items():
        reverse_reactant_map[value] = key
    for key, value in product_map.items():
        reverse_product_map[value] = key
    return reactant_map, product_map, reverse_reactant_map, reverse_product_map


def find_center(smarts, automapped_smarts, maps):
    def find_max_index_in_smarts(automapped_smarts):
        pattern_str = r"\[[A-Z0-9]+:(\d+)\]"
        pattern = re.compile(pattern_str)
        iter = re.finditer(pattern, automapped_smarts)
        max_index = -1
        for matcher in iter:
            try:
                index = int(matcher.group(1))
                if index > max_index:
                    max_index = index
            except Exception as e:
                print(e)
        if max_index >= 0:
            return max_index
        else:
            raise Exception("No effective index")

    def compare_both_side(index, maps, ractant_matrices, product_matrices):
        """

        :param index: 要比较的索引值
        :param maps: 4个map分别是reactant_map, product_map, reverse_reactant_map, reverse_product_map
        :param ractant_matrices: 反应物的matrix，由get_adjacent_matrix方法读取
        :param product_matrices: 产物的matrix，由get_adjacent_matrix方法读取
        """
        assert len(maps) == 4
        reactant_map, product_map, reverse_reactant_map, reverse_product_map = maps[0], maps[1], maps[2], maps[3]
        r_id, r_index = reactant_map.get(str(index)).split("_")
        r_id, r_index = int(r_id), int(r_index)
        p_id, p_index = product_map.get(str(index)).split("_")
        p_id, p_index = int(p_id), int(p_index)
        r_M = ractant_matrices[r_id]
        p_M = product_matrices[p_id]
        r_row = r_M[r_index, :].tolist()
        p_row = p_M[p_index, :].tolist()
        r_map = dict()
        p_map = dict()
        for i in range(len(r_row)):
            if r_row[i] >= 1:
                record_index = "{:d}_{:d}".format(r_id, i)
                automapped_index = reverse_reactant_map.get(record_index)
                r_map[automapped_index] = r_row[i]
        for i in range(len(p_row)):
            if p_row[i] >= 1:
                record_index = "{:d}_{:d}".format(p_id, i)
                automapped_index = reverse_product_map.get(record_index)
                p_map[automapped_index] = p_row[i]
        return r_map != p_map

    # def merge_centers(centers, map, smiles, dist):

    def classify_and_merge_indices(indices, matrices):
        result_map = dict()
        for index in indices:
            atom_id, atom_pos = index.split("_")
            if result_map.get(str(atom_id)) is None:
                result_map[str(atom_id)] = dict(changed_pos=[index])
            else:
                result_map[str(atom_id)]["changed_pos"].append(index)
        for key, values in result_map.items():
            M = matrices[int(key)]
            point_set = set()
            for index in values["changed_pos"]:
                sub_matrix, points = get_neighbor_graph(M, int(index.split("_")[1]))
                for point in points:
                    point_set.add(point)
            result_map[key]["merged_center"] = point_set
        return result_map

    parts = smarts.split(">>")
    lhs = parts[0]
    rhs = parts[1]
    reactants = lhs.split(".")
    products = rhs.split(".")
    reactant_matrices = []
    product_matrices = []
    for reactant in reactants:
        M, mol = get_adjacent_matrix(reactant)
        reactant_matrices.append(M)
    for product in products:
        M, mol = get_adjacent_matrix(product)
        product_matrices.append(M)

    max_index = find_max_index_in_smarts(automapped_smarts)
    centers = []
    for index in range(1, max_index + 1):
        if compare_both_side(index, maps, reactant_matrices, product_matrices):
            centers.append(index)
    reactant_indices = []
    product_indices = []
    for center in centers:
        reactant_indices.append(maps[0][str(center)])
        product_indices.append(maps[1][str(center)])
    reactant_center_map = classify_and_merge_indices(reactant_indices, reactant_matrices)
    product_center_map = classify_and_merge_indices(product_indices, product_matrices)
    return reactant_center_map, product_center_map


def get_rule(smarts, mapped_smarts, reactant_center_map=None, product_center_map=None):
    reactant_map, product_map, reverse_reactant_map, reverse_product_map = get_index_map_from_smarts(mapped_smarts)
    if reactant_center_map is None or product_center_map is None:
        reactant_center_map, product_center_map = find_center(smarts, mapped_smarts,
                                                              [reactant_map, product_map, reverse_reactant_map,
                                                               reverse_product_map])
    for key in reactant_center_map.keys():
        print(key)
        print(reactant_center_map[key]["changed_pos"])
        print(reactant_center_map[key]["merged_center"])
    print("#" * 20)
    for key in product_center_map.keys():
        print(key)
        print(product_center_map[key]["changed_pos"])
        print(product_center_map[key]["merged_center"])
    parts = smarts.split(">>")
    lhs = parts[0]
    rhs = parts[1]
    reactants = lhs.split(".")
    products = rhs.split(".")
    reac_num = len(reactants)
    prod_num = len(products)
    '''temporary'''
    assert reac_num == 2
    assert prod_num == 1
    '''temporary'''
    assert len(reactant_center_map) == reac_num
    assert len(product_center_map) == prod_num
    sub_smiles = []
    reactant_sub_smiles_list = []
    product_sub_smiles_list = []
    # reactant parts
    for i, molecule in enumerate(reactants):
        points = reactant_center_map.get(str(i)).get("merged_center")
        reactant_sub_smiles = get_sub_smiles(molecule, points)
        reactant_sub_smiles_list.append(reactant_sub_smiles)
    # produtc parts
    for i, molecule in enumerate(products):
        points = product_center_map.get(str(i)).get("merged_center")
        product_sub_smiles = get_sub_smiles(molecule, points)
        product_sub_smiles_list.append(product_sub_smiles)
    rule = ".".join(reactant_sub_smiles_list) + ">>" + ".".join(product_sub_smiles_list)
    return rule


if __name__ == "__main__":
    # # smiles = "C1(C)(C)c2c(cccc2)N(c3ccc(N(c4ccc(N5c6c(cccc6)C(C)(C)c(cccc7)c57)cc4)c8ccc([nH]c(ccc(c9ccccc9)c%10)c%10%11)c%11c8)cc3)c(cccc%12)c1%12"
    # # smiles = "c12c(ccc(Br)c1)C=Cc(cccc3)c3C2=O"
    # smiles = "C([Mg]Br)CCC.C([Mg]Br)CCC"
    # # smiles = "C(CCCCC)CCCCN"
    # smarts = "[CH2:1]([CH2:3][CH2:5][CH2:7][CH2:9][CH3:11])[CH2:2][CH2:4][CH2:6][CH2:8][NH2:10].[CH2:12]([CH2:14][CH2:16][CH2:18][CH2:20][CH3:22])[CH2:13][CH2:15][CH2:17][CH2:19][NH2:21]>>C(=O)([NH:21][CH2:19][CH2:17][CH2:15][CH2:13][CH2:12][CH2:14][CH2:16][CH2:18][CH2:20][CH3:22])[NH:10][CH2:8][CH2:6][CH2:4][CH2:2][CH2:1][CH2:3][CH2:5][CH2:7][CH2:9][CH3:11]"
    # # C([Mg]Br)CCC.C([Mg]Br)CCC>>C(/C(C)=O)=C(/O)\CCCC
    # smarts1 = "[CH2:1]([Mg]Br)[CH2:2][CH2:4][CH3:6].[CH2:7]([Mg]Br)[CH2:8][CH2:10][CH3:12]>>[CH:2](/[C:4]([CH3:6])=O)=[C:1](/O)\[CH2:7][CH2:8][CH2:10][CH3:12]"
    # get_sub_smiles(smiles, 0)
    # result = get_index_map_from_smarts(smarts1)
    # print(result)
    smarts = "C(CCCCC)CCCCN.C(CCCCC)CCCCN>>C(=O)(NCCCCCCCCCC)NCCCCCCCCCC"
    mapped_smarts = "[CH2:1]([CH2:3][CH2:5][CH2:7][CH2:9][CH3:11])[CH2:2][CH2:4][CH2:6][CH2:8][NH2:10].[CH2:12]([CH2:14][CH2:16][CH2:18][CH2:20][CH3:22])[CH2:13][CH2:15][CH2:17][CH2:19][NH2:21]>>C(=O)([NH:21][CH2:19][CH2:17][CH2:15][CH2:13][CH2:12][CH2:14][CH2:16][CH2:18][CH2:20][CH3:22])[NH:10][CH2:8][CH2:6][CH2:4][CH2:2][CH2:1][CH2:3][CH2:5][CH2:7][CH2:9][CH3:11]"
    s = "C(C#C)N(C)C.C(=O)=O>>N(C)(C)CC#CC(=O)O"
    m = "[CH2:1]([C:3]#[CH:6])[N:2]([CH3:5])[CH3:4].[C:7](=[O:9])=[O:8]>>[N:2]([CH3:5])([CH3:4])[CH2:1][C:3]#[C:6][C:7](=[O:9])[OH:8]"
    print("规则\n", get_rule(smarts, mapped_smarts))
