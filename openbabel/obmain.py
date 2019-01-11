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
import json
import os
import re

import mysql.connector
import numpy as np
import openbabel as ob
from rdkit import Chem
from rdkit.Chem import AllChem

from dataop import read_all_data


class SkipException(Exception):
    pass


class NotFormalException(Exception):
    pass


class NotFormalConvertion(Exception):
    pass


def get_adjacent_matrix(smiles):
    """
    用openbabel读smiles
    :param smiles:smiles字符串
    :return:
        M是邻接矩阵
        mol是读smiles之后的内部对象
    """
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
    """
    功能同 @{function}get_adjacent_matrix，中间通过sdf文件格式转化
    :param smiles:
    :return:
    """
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
    return M, mol


def get_neighbor_graph(matrix, index, dist=3):
    """
    根据邻接矩阵matrix和原子索引index获取dist距离的中心
    :param matrix: 邻接矩阵
    :param index: 原子索引
    :param dist: 半径或者说步长
    :return: 中心邻接矩阵和其中对应的点在之前的矩阵中的序号
    """
    assert type(matrix) == np.ndarray
    if index > matrix.shape[0] or index < 0:
        raise ValueError(index)
    neighbors = dict()
    neighbors["neighbor_1"] = set()
    neighbors["neighbor_1"].add(index)
    index_list = np.array(range(matrix.shape[0]))
    target_points = set()
    target_points.add(index)
    for i in range(2, dist + 1):
        prev_key = "neighbor_{:d}".format(i - 1)
        cur_key = "neighbor_{:d}".format(i)
        neighbors[cur_key] = set()
        prev_points = neighbors[prev_key]
        for point in prev_points:
            for adj_point in index_list[matrix[point, :] > 0].tolist():
                if adj_point not in target_points:
                    neighbors[cur_key].add(adj_point)
                    target_points.add(adj_point)
    target_points = list(target_points)
    target_points.sort()
    return matrix[target_points, :][:, target_points], target_points


def get_neighbor_graph_exclude_aromatic_C(matrix, mol, index):
    """
    功能类同@{function}get_neighbor_graph，但是会排除芳香的C
    """
    assert type(mol) == ob.OBMol
    assert type(matrix) == np.ndarray
    if index > matrix.shape[0] or index < 0:
        raise ValueError(index)
    neighbors = dict()
    neighbors["neighbor_1"] = set()
    neighbors["neighbor_1"].add(index)
    index_list = np.array(range(matrix.shape[0]))
    target_points = set()
    target_points.add(index)
    # 1步中心已经设置好，从2开始
    steps = 2
    while len(target_points) < matrix.shape[0] and steps < 1000:
        prev_key = "neighbor_{:d}".format(steps - 1)
        cur_key = "neighbor_{:d}".format(steps)
        neighbors[cur_key] = set()
        prev_points = neighbors[prev_key]
        for point in prev_points:
            for new_point in index_list[matrix[point, :] > 0].tolist():
                if new_point not in target_points:
                    # 过滤非芳香性的碳原子C
                    atom = mol.GetAtomById(new_point)
                    if atom.IsCarbon() and not atom.IsAromatic():
                        continue
                    neighbors[cur_key].add(new_point)
                    target_points.add(point)
        steps = steps + 1
    if steps == 1000:
        raise ValueError("Too many extend steps")
    target_points = list(target_points)
    target_points.sort()
    return matrix[target_points, :][:, target_points], target_points


def get_sub_smiles(smiles, points=None, pos=None, dist=3, Exclude_Aromatic_C=False):
    # print("dist=", dist)
    M, mol = get_adjacent_matrix(smiles)
    # print("分子内部顺序")
    # for atom in ob.OBMolAtomIter(mol):
    #     print(atom.GetId(), atom.GetIdx(), atom.GetIndex(), atom.GetAtomicNum(), atom.GetType())
    if points is None and pos is not None:
        if Exclude_Aromatic_C:
            sub_M, points = get_neighbor_graph_exclude_aromatic_C(M, mol, pos)
        else:
            sub_M, points = get_neighbor_graph(M, pos, dist)

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
    # print("中心分子内部顺序")
    # for atom in ob.OBMolAtomIter(new_mol):
    #     print(atom.GetId(), atom.GetIdx(), atom.GetIndex(), atom.GetAtomicNum(), atom.GetType())
    OBConversion.SetOutFormat("smi")
    new_mol_smiles = OBConversion.WriteString(new_mol).strip()
    # print("中心smiles")
    # print(new_mol_smiles)
    # print("#" * 80)
    return new_mol_smiles


def get_index_map_from_smarts(mapped_smarts):
    assert len(mapped_smarts.split(">>")) == 2
    error_msg = []
    # left hand side   right hand side
    parts = mapped_smarts.split(">>")
    lhs = parts[0]
    rhs = parts[1]
    reactants = lhs.split(".")
    products = rhs.split(".")
    pattern_str = r"\[[A-Za-z0-9]+:(\d+)\]|(Mg|Br|[A-Z]+[a-z]*)"
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
                error_msg.append("{:d}号位原子{:s}没有匹配，在这个smarts中：{:s}".format(count, matcher.group(2), reactant))

            count = count + 1
    for i, product in enumerate(products):
        count = 0
        iter = re.finditer(pattern, product)
        for matcher in iter:
            automap_index = matcher.group(1)
            if automap_index:
                product_map[str(automap_index)] = get_record_index(i, count)
            else:
                error_msg.append("{:d}号位原子{:s}没有匹配，在这个smarts中：{:s}".format(count, matcher.group(2), product))
            count = count + 1
    reverse_reactant_map = dict()
    reverse_product_map = dict()
    for key, value in reactant_map.items():
        reverse_reactant_map[value] = key
    for key, value in product_map.items():
        reverse_product_map[value] = key
    return reactant_map, product_map, reverse_reactant_map, reverse_product_map, error_msg


def find_center(smarts, automapped_smarts, maps, dist=3, Exclude_Aromatic_C=False):
    def find_max_index_in_smarts(automapped_smarts):
        """
        找到automap里最大的index值
        :param automapped_smarts:
        :return:
        """
        pattern_str = r"\[[A-Za-z0-9]+:(\d+)\]"
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
            raise ValueError("No effective index")

    def compare_both_side(index, maps, reactant_matrices, product_matrices):
        """
        比较两端同一个index的环境是否发生变化
        :param index: 要比较的索引值(automap序号)
        :param maps: 4个map分别是reactant_map, product_map, reverse_reactant_map, reverse_product_map
        :param ractant_matrices: 反应物的matrix，由get_adjacent_matrix方法读取
        :param product_matrices: 产物的matrix，由get_adjacent_matrix方法读取
        """
        assert len(maps) == 4
        reactant_map, product_map, reverse_reactant_map, reverse_product_map = maps[0], maps[1], maps[2], maps[3]
        # 通过index找到反应物中对应的{反应物序号}和{原子序号}
        r_id, r_index = reactant_map.get(str(index)).split("_")
        r_id, r_index = int(r_id), int(r_index)
        # 通过index找到产物中对应的{反应物序号}和{原子序号}
        p_id, p_index = product_map.get(str(index)).split("_")
        p_id, p_index = int(p_id), int(p_index)
        # 找到对应矩阵
        r_M = reactant_matrices[r_id]
        p_M = product_matrices[p_id]
        # 找到矩阵中原子对应的行（有哪些邻接点）
        r_row = r_M[r_index, :].tolist()
        p_row = p_M[p_index, :].tolist()

        r_map = dict()
        r_map["UnknownAtoms"] = set()
        p_map = dict()
        p_map["UnknownAtoms"] = set()
        for i in range(len(r_row)):
            # 度数大于1即有键
            if r_row[i] >= 1:
                record_index = "{:d}_{:d}".format(r_id, i)
                automapped_index = reverse_reactant_map.get(record_index)
                if automapped_index:
                    r_map[automapped_index] = r_row[i]
                else:
                    r_map["UnknownAtoms"].add(automapped_index)

        for i in range(len(p_row)):
            if p_row[i] >= 1:
                record_index = "{:d}_{:d}".format(p_id, i)
                automapped_index = reverse_product_map.get(record_index)
                if automapped_index:
                    p_map[automapped_index] = p_row[i]
                else:
                    p_map["UnknownAtoms"].add(automapped_index)
        return r_map != p_map

    def classify_and_merge_indices(indices, matrices, mols, dist=3, Exclude_Aromatic_C=False):
        result_map = dict()
        for index in indices:
            atom_id, atom_pos = index.split("_")
            if result_map.get(str(atom_id)) is None:
                result_map[str(atom_id)] = dict(changed_pos=[int(atom_pos)])
            else:
                result_map[str(atom_id)]["changed_pos"].append(int(atom_pos))
        for key, values in result_map.items():
            M = matrices[int(key)]
            mol = mols[int(key)]
            point_set = set()
            for index in values["changed_pos"]:
                if Exclude_Aromatic_C:
                    sub_matrix, points = get_neighbor_graph_exclude_aromatic_C(M, mol, int(index))
                else:
                    # 每一个改变的点，都扩展dist步，将它们加到point_set里
                    sub_matrix, points = get_neighbor_graph(M, int(index), dist)
                for point in points:
                    point_set.add(point)
            result_map[key]["merged_center"] = point_set
        return result_map

    def connect_shortest_path(points, matrix):
        """
        给一个矩阵和两个点，找它们之间的最短路径
        :param points:
        :param matrix:
        :return:
        """
        if len(points) != 2:
            return points

        def convert_dist_matrix(M):
            assert isinstance(M, np.ndarray)
            M = np.int32(M)
            M[M > 0] = 1
            return M

        i, j = points[0], points[1]
        M = convert_dist_matrix(matrix)
        point_num = M.shape[0]
        dist = np.squeeze(M[i])
        pre_node_list = [-1] * point_num
        visited = set()
        visited.add(i)
        unreached = set(range(0, point_num))
        unreached.remove(i)
        while j not in visited:
            pre_node = -1
            min_dist = 10e6
            k = -1
            for v in visited:
                for neighbor_v in np.array(range(0, point_num))[M[v] > 0]:
                    if neighbor_v not in visited:
                        if M[v, neighbor_v] < min_dist:
                            min_dist = M[v, neighbor_v]
                            pre_node = v
                            k = neighbor_v
            if k == -1:
                break
            visited.add(k)
            unreached.remove(k)
            for index in range(0, point_num):
                if M[k, index] > 0 and dist[k] + M[k, index] < dist[index]:
                    dist[index] = dist[k] + M[k, index]
            pre_node_list[k] = pre_node
        if pre_node_list[j] != -1:
            path_nodes = []
            temp = j
            while pre_node_list[temp] != -1:
                path_nodes.append(temp)
                temp = pre_node_list[temp]
            path_nodes.append(i)
            return path_nodes
        else:
            return points

    error_msg = []

    parts = smarts.split(">>")
    lhs = parts[0]
    rhs = parts[1]
    reactants = lhs.split(".")
    products = rhs.split(".")
    reactant_matrices = []
    product_matrices = []
    reactant_obmols = []
    product_obmols = []
    for reactant in reactants:
        M, mol = get_adjacent_matrix(reactant)
        reactant_matrices.append(M)
        reactant_obmols.append(mol)
    for product in products:
        M, mol = get_adjacent_matrix(product)
        product_matrices.append(M)
        product_obmols.append(mol)

    max_index = find_max_index_in_smarts(automapped_smarts)
    centers = []
    for index in range(1, max_index + 1):
        try:
            if compare_both_side(index, maps, reactant_matrices, product_matrices):
                centers.append(index)
        except AttributeError as e:
            msg = "编号为{:d}的原子在这个smarts中没出现：{:s}".format(index, automapped_smarts)
            error_msg.append(msg)
    reactant_indices = []
    product_indices = []
    for center in centers:
        reactant_indices.append(maps[0][str(center)])
        product_indices.append(maps[1][str(center)])

    """同一个分子的信息集中到一起，放到一个map里"""
    reactant_center_map = classify_and_merge_indices(reactant_indices, reactant_matrices, reactant_obmols, dist,
                                                     Exclude_Aromatic_C=Exclude_Aromatic_C)
    product_center_map = classify_and_merge_indices(product_indices, product_matrices, reactant_obmols, dist,
                                                    Exclude_Aromatic_C=Exclude_Aromatic_C)

    """添加最短路径"""
    for key in reactant_center_map.keys():
        points = reactant_center_map[key]["changed_pos"]
        if len(points) == 2:
            temp_set = reactant_center_map[key]["merged_center"]
            for path_point in connect_shortest_path(points, reactant_matrices[int(key)]):
                temp_set.add(path_point)
            reactant_center_map[key]["merged_center"] = temp_set
    for key in product_center_map.keys():
        points = product_center_map[key]["changed_pos"]
        if len(points) == 2:
            temp_set = product_center_map[key]["merged_center"]
            for path_point in connect_shortest_path(points, product_matrices[int(key)]):
                temp_set.add(path_point)
            product_center_map[key]["merged_center"] = temp_set

    return reactant_center_map, product_center_map, error_msg


def check_rule(rule):
    """
    检查规则反应物，产物数量
    :param rule:
    :return:
    """
    try:
        parts = rule.split(">>")
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
    except AssertionError:
        return False
    return True


def get_rule(smarts, mapped_smarts, reactant_center_map=None, product_center_map=None, dist=3,
             Exclude_Aromatic_C=False):
    """
    获取规则
    :param smarts:
    :param mapped_smarts:
    :param reactant_center_map:
    :param product_center_map:
    :param dist:
    :return:
    """
    map_error_msg = []
    find_center_error_msg = []
    if reactant_center_map is None or product_center_map is None:
        reactant_map, product_map, reverse_reactant_map, reverse_product_map, map_error_msg = get_index_map_from_smarts(mapped_smarts)
        reactant_center_map, product_center_map, find_center_error_msg = find_center(smarts, mapped_smarts, [reactant_map, product_map, reverse_reactant_map, reverse_product_map], dist,
                                                                                     Exclude_Aromatic_C=Exclude_Aromatic_C)

    # for key in reactant_center_map.keys():
    #     print(key)
    #     print(reactant_center_map[key]["changed_pos"])
    #     print(reactant_center_map[key]["merged_center"])
    # print("#" * 20)
    # for key in product_center_map.keys():
    #     print(key)
    #     print(product_center_map[key]["changed_pos"])
    #     print(product_center_map[key]["merged_center"])

    # 检查原本的规则是不是2to1
    if not check_rule(smarts):
        raise NotFormalException

    parts = smarts.split(">>")
    lhs = parts[0]
    rhs = parts[1]
    reactants = lhs.split(".")
    products = rhs.split(".")
    reac_num = len(reactants)
    prod_num = len(products)
    assert len(reactant_center_map) == reac_num
    assert len(product_center_map) == prod_num
    reactant_sub_smiles_list = []
    product_sub_smiles_list = []

    # reactant parts
    for i, molecule in enumerate(reactants):
        points = reactant_center_map.get(str(i)).get("merged_center")
        reactant_sub_smiles = get_sub_smiles(molecule, points=points, dist=1, Exclude_Aromatic_C=Exclude_Aromatic_C)
        reactant_sub_smiles_list.append(reactant_sub_smiles)
    # produtc parts
    for i, molecule in enumerate(products):
        points = product_center_map.get(str(i)).get("merged_center")
        product_sub_smiles = get_sub_smiles(molecule, points=points, dist=1, Exclude_Aromatic_C=Exclude_Aromatic_C)
        product_sub_smiles_list.append(product_sub_smiles)
    rule = ".".join(reactant_sub_smiles_list) + ">>" + ".".join(product_sub_smiles_list)

    # 检查提取后是不是2to1
    if not check_rule(rule):
        raise NotFormalConvertion
    all_error_msg = map_error_msg + find_center_error_msg
    # for msg in all_error_msg:
    #     print(msg, file=sys.stderr)
    return rule, reactant_center_map, product_center_map, all_error_msg


def get_dict(smarts, mapped_smarts):
    """
    获取打印信息的字典
    :param smarts:smiles格式的反应
    :param mapped_smarts:map之后的反应
    :return:
    """
    dictionary = dict()
    dictionary["Original_Smiles"] = smarts
    dictionary["Smiles"] = mapped_smarts
    dictionary["Distance"] = []
    reactant_map, product_map, reverse_reactant_map, reverse_product_map, error_msg = get_index_map_from_smarts(
        mapped_smarts)
    for dist in (1, 2, 3):
        rule, r, p, error_msg = get_rule(smarts, mapped_smarts, dist=dist)
        if not check_rule(rule):
            return None
        temp_dict = dict()
        temp_dict["Centers"] = dict(Reactants=dict(), Products=dict())
        temp_dict["Centers2"] = dict(Reactants=dict(), Products=dict())
        temp_dict["Error_message"] = error_msg
        for key in r.keys():
            temp_dict["Centers"]["Reactants"][key] = str(r[key]["merged_center"])
            indices = []
            for pos in r[key]["merged_center"]:
                new_key = "{:s}_{:s}".format(key, str(pos))
                index = reverse_reactant_map.get(new_key)
                indices.append(index)
            temp_dict["Centers2"]["Reactants"][key] = str(indices)
        for key in p.keys():
            temp_dict["Centers"]["Products"][key] = str(p[key]["merged_center"])
            indices = []
            for pos in p[key]["merged_center"]:
                new_key = "{:s}_{:s}".format(key, str(pos))
                index = reverse_product_map.get(new_key)
                indices.append(index)
            temp_dict["Centers2"]["Products"][key] = str(indices)
        temp_dict["Rule"] = rule
        dictionary["Distance"].append({str(dist): temp_dict})
        if dist == 3:
            dictionary["Changed_atoms"] = dict(Reactants=dict(), Products=dict())
            dictionary["Changed_atoms2"] = dict(Reactants=dict(), Products=dict())
            for key in r.keys():
                dictionary["Changed_atoms"]["Reactants"][key] = str(r[key]["changed_pos"])
                indices = []
                for pos in r[key]["changed_pos"]:
                    new_key = "{:s}_{:s}".format(key, str(pos))
                    index = reverse_reactant_map.get(new_key)
                    indices.append(index)
                dictionary["Changed_atoms2"]["Reactants"][key] = str(indices)
            for key in p.keys():
                dictionary["Changed_atoms"]["Products"][key] = str(p[key]["changed_pos"])
                indices = []
                for pos in p[key]["changed_pos"]:
                    new_key = "{:s}_{:s}".format(key, str(pos))
                    index = reverse_product_map.get(new_key)
                    indices.append(index)
                dictionary["Changed_atoms2"]["Products"][key] = str(indices)
    return dictionary


def smiles2fp(smiles):
    fp = None
    try:
        mol = Chem.MolFromSmiles(smiles)
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048).ToBitString()
    except Exception as e:
        print(e)
    finally:
        return fp


def main1():
    data = read_all_data()
    count = 1
    num = 1
    AE, TE, VE, OE = 0, 0, 0, 0
    for l in data:
        error_msg = ""
        oe = ""
        result = ""
        i = 1
        for item in l:
            try:
                osmi, smi = item[0], item[1]
                d = get_dict(osmi, smi)
                if d:
                    result += json.dumps(d, indent=4, sort_keys=True, ensure_ascii=False) + "\n"
            except AssertionError as e:
                AE = AE + 1
                error_msg = error_msg + (str(item) + "\n" + str(e) + "\n")
            except TypeError as e:
                TE = TE + 1
                error_msg = error_msg + (str(item) + "\n" + str(e) + "\n")
            except ValueError as e:
                VE = VE + 1
                error_msg = error_msg + (str(item) + "\n" + str(e) + "\n")
            except Exception as e:
                OE = OE + 1
                oe = oe + (str(item) + "\n" + str(e) + "\n")
            finally:
                print("{:d}_{:d}".format(count, i))
                i = i + 1
                num = num + 1
        path = os.path.abspath(os.path.join("..", "rule", "result_{:d}.txt".format(count)))
        path2 = os.path.abspath(os.path.join("..", "rule", "error_{:d}.txt".format(count)))
        if not os.path.isdir(os.path.dirname(path)):
            os.mkdir(os.path.dirname(path))
        f = open(path, "w")
        f.write(result)
        f.close()
        f = open(path2, "w")
        f.write(error_msg)
        f.close()
        f = open(path2, "a+")
        f.write(oe)
        f.close()
        count = count + 1
    print("Num:{:d},AssertionError:{:d},TypeError:{:d},ValueError:{:d},OtherError:{:d}".format(num, AE, TE, VE, OE))
    # s1 = "C(CCCCC)CCCCN.C(CCCCC)CCCCN>>C(=O)(NCCCCCCCCCC)NCCCCCCCCCC"
    # m1 = "[CH2:1]([CH2:3][CH2:5][CH2:7][CH2:9][CH3:11])[CH2:2][CH2:4][CH2:6][CH2:8][NH2:10].[CH2:12]([CH2:14][CH2:16][CH2:18][CH2:20][CH3:22])[CH2:13][CH2:15][CH2:17][CH2:19][NH2:21]>>C(=O)([NH:21][CH2:19][CH2:17][CH2:15][CH2:13][CH2:12][CH2:14][CH2:16][CH2:18][CH2:20][CH3:22])[NH:10][CH2:8][CH2:6][CH2:4][CH2:2][CH2:1][CH2:3][CH2:5][CH2:7][CH2:9][CH3:11]"
    # s2 = "C(C#C)N(C)C.C(=O)=O>>N(C)(C)CC#CC(=O)O"
    # m2 = "[CH2:1]([C:3]#[CH:6])[N:2]([CH3:5])[CH3:4].[C:7](=[O:9])=[O:8]>>[N:2]([CH3:5])([CH3:4])[CH2:1][C:3]#[C:6][C:7](=[O:9])[OH:8]"
    # s3 = "c1(ccc(C)c(N)c1)C(=O)OC.c1(C(=O)Cl)cc(OC)ccc1Br>>c1(cc(OC)ccc1Br)C(=O)Nc2cc(C(=O)OC)ccc2C"
    # m3 = "[c:1]1([cH:4][cH:8][c:9]([CH3:12])[c:5]([NH2:10])[cH:2]1)[C:3](=[O:7])[O:6][CH3:11].[c:13]1([C:16](=[O:21])Cl)[cH:15][c:19]([O:23][CH3:24])[cH:22][cH:17][c:14]1[Br:18]>>[c:13]1([cH:15][c:19]([O:23][CH3:24])[cH:22][cH:17][c:14]1[Br:18])[C:16](=[O:21])[NH:10][c:5]2[cH:2][c:1]([C:3](=[O:7])[O:6][CH3:11])[cH:4][cH:8][c:9]2[CH3:12]"
    # s4 = "C(C[Hg]Br)CCCC.C(C[Hg]Br)CCCC>>C(=O)(CCCCCC)CCCCCC"
    # m4 = "[CH2:1]([CH2:3][Hg]Br)[CH2:2][CH2:4][CH2:6][CH3:8].[CH2:9]([CH2:11][Hg]Br)[CH2:10][CH2:12][CH2:14][CH3:16]>>C(=O)([CH2:11][CH2:9][CH2:10][CH2:12][CH2:14][CH3:16])[CH2:3][CH2:1][CH2:2][CH2:4][CH2:6][CH3:8]"
    # s5 = "C(=O)(Cl)CCCCCCCCC.C(=O)(Cl)CCCCCCCCC>>C(O)(COC(=O)CCCCCCCCC)COC(=O)CCCCCCCCC"
    # m5 = "[C:1](=[O:4])(Cl)[CH2:2][CH2:5][CH2:6][CH2:7][CH2:8][CH2:9][CH2:10][CH2:11][CH3:12].[C:13](=[O:16])(Cl)[CH2:14][CH2:17][CH2:18][CH2:19][CH2:20][CH2:21][CH2:22][CH2:23][CH3:24]>>C(O)(CO[C:13](=[O:16])[CH2:14][CH2:17][CH2:18][CH2:19][CH2:20][CH2:21][CH2:22][CH2:23][CH3:24])CO[C:1](=[O:4])[CH2:2][CH2:5][CH2:6][CH2:7][CH2:8][CH2:9][CH2:10][CH2:11][CH3:12]"
    # s6 = "C(C)O.C(CC)(CCC(=O)OCC)N(=O)=O>>N1(CC)CCCC1CC"
    # m6 = "[CH2:1]([CH3:3])O.[CH:4]([CH2:7][CH3:11])([CH2:6][CH2:10][C:12](=O)OCC)[N:5](=O)=O>>[N:5]1([CH2:1][CH3:3])[CH2:12][CH2:10][CH2:6][CH:4]1[CH2:7][CH3:11]"
    # s7 = "C(C)(C)(N)CC.Cl[H]>>C(C)(C)(CC)N=C=S"
    # m7 = "[C:1]([CH3:5])([CH3:4])([NH2:3])[CH2:2][CH3:6].Cl[H]>>[C:1]([CH3:5])([CH3:4])([CH2:2][CH3:6])[N:3]=C=S"
    # s8 = "C(CBr)CCCl.C(CBr)CCCl>>C(CCCCCl)#CCCCCCl"
    # m8 = "[CH2:1]([CH2:3]Br)[CH2:2][CH2:4][Cl:6].[CH2:7]([CH2:9]Br)[CH2:8][CH2:10][Cl:12]>>C([CH2:3][CH2:1][CH2:2][CH2:4][Cl:6])#C[CH2:9][CH2:7][CH2:8][CH2:10][Cl:12]"
    # s9 = "C(N)C=C.C(N)C=C>>C(=O)(NCC=C)OCCCCl"
    # m9 = "[CH2:1]([NH2:3])[CH:2]=[CH2:4].[CH2:5](N)[CH:6]=[CH2:8]>>C(=O)([NH:3][CH2:1][CH:2]=[CH2:4])O[CH2:5][CH2:6][CH2:8]Cl"
    # s10 = "C(CCC)CCI.C(CCC)CCI>>C(C#N)(CCCCCC)(NC(C)=O)C(=O)OCC"
    # m10 = "[CH2:1]([CH2:3][CH2:5][CH3:7])[CH2:2][CH2:4]I.[CH2:8]([CH2:10][CH2:12][CH3:14])[CH2:9][CH2:15]I>>[C:8]([C:10]#N)([CH2:4][CH2:2][CH2:1][CH2:3][CH2:5][CH3:7])(NC([CH3:15])=O)[C:9](=O)O[CH2:12][CH3:14]"
    # s11 = "C(CCC)CCN.N(CC)CC>>C(=O)(NCCCCCC)N(CC)CC"
    # m11 = "[CH2:1]([CH2:3][CH2:5][CH3:7])[CH2:2][CH2:4][NH2:6].[NH:8]([CH2:10][CH3:12])[CH2:9][CH3:11]>>C(=O)([NH:6][CH2:4][CH2:2][CH2:1][CH2:3][CH2:5][CH3:7])[N:8]([CH2:10][CH3:12])[CH2:9][CH3:11]"
    # s12 = "C(C)(C)CNCC(C)C.C(CCC)CCN>>N(CC(C)C)(CC(C)C)C(=O)NCCCCCC"
    # m12 = "[CH:1]([CH3:4])([CH3:3])[CH2:2][NH:5][CH2:6][CH:7]([CH3:9])[CH3:8].[CH2:10]([CH2:12][CH2:14][CH3:16])[CH2:11][CH2:13][NH2:15]>>[N:5]([CH2:6][CH:7]([CH3:9])[CH3:8])([CH2:2][CH:1]([CH3:4])[CH3:3])C(=O)[NH:15][CH2:13][CH2:11][CH2:10][CH2:12][CH2:14][CH3:16]"
    # s13 = "C(=O)(O)CCl.C(=O)(O)CCl>>C(O)(CCl)COC(=O)CCl"
    # m13 = "[C:1](=[O:4])([OH:3])[CH2:2][Cl:5].[C:6](=O)([OH:8])[CH2:7][Cl:10]>>[CH:6]([OH:8])([CH2:7][Cl:10])C[O:3][C:1](=[O:4])[CH2:2][Cl:5]"
    # s14 = "C(C)(=O)Cl.C(C)(=O)Cl.C(C)(=O)Cl>>C(C#N)(C(C)=O)C(=O)OC"
    # m14 = "[C:1]([CH3:4])(=[O:3])Cl.[C:5]([CH3:8])(=[O:7])Cl.[C:9]([CH3:12])(=[O:11])Cl>>[CH:8]([C:12]#N)([C:1]([CH3:4])=[O:3])[C:5](=[O:7])[O:11][CH3:9]"
    # s15 = "N(C)C.N(C)C.N(C)C>>C(=N)(NC(=N)N(C)C)N(C)C"
    # m15 = "[NH:1]([CH3:3])[CH3:2].[NH:4]([CH3:6])[CH3:5].[NH:7]([CH3:9])[CH3:8]>>[C:8](=N)([NH:7][C:9](=N)[N:4]([CH3:6])[CH3:5])[N:1]([CH3:3])[CH3:2]"
    # s16 = "C(CN)CCC.C(CN)CCC>>C(O)(CC(=O)NCCCCC)C(=O)NCCCCC"
    # m16 = "[CH2:1]([CH2:3][NH2:5])[CH2:2][CH2:4][CH3:6].[CH2:7]([CH2:9][NH2:11])[CH2:8][CH2:10][CH3:12]>>C(O)(CC(=O)[NH:11][CH2:9][CH2:7][CH2:8][CH2:10][CH3:12])C(=O)[NH:5][CH2:3][CH2:1][CH2:2][CH2:4][CH3:6]"
    # s17 = "C(=O)C=C.C(=O)C=C>>C(O)(C=C)CCCOC"
    # m17 = "[CH:1](=[O:3])[CH:2]=[CH2:4].[CH:5](=[O:7])[CH:6]=[CH2:8]>>[CH:1]([OH:3])([CH:2]=[CH2:4])[CH2:8][CH2:6][CH2:5][O:7]C"
    # s18 = "N(CCCC)CCCC.C#C>>N(CCCC)(CCCC)C(=O)C=C"
    # m18 = "[NH:1]([CH2:3][CH2:5][CH2:7][CH3:9])[CH2:2][CH2:4][CH2:6][CH3:8].[CH:10]#[CH:11]>>[N:1]([CH2:3][CH2:5][CH2:7][CH3:9])([CH2:2][CH2:4][CH2:6][CH3:8])C(=O)[CH:10]=[CH2:11]"
    # s19 = "C(CCC)CCBr.C(CCC)CCBr>>S(=O)(=O)(CCCCCC)CCCCCC"
    # m19 = "[CH2:1]([CH2:3][CH2:5][CH3:7])[CH2:2][CH2:4]Br.[CH2:8]([CH2:10][CH2:12][CH3:14])[CH2:9][CH2:11]Br>>S(=O)(=O)([CH2:11][CH2:9][CH2:8][CH2:10][CH2:12][CH3:14])[CH2:4][CH2:2][CH2:1][CH2:3][CH2:5][CH3:7]"
    # s20 = "C(CCCCCCCCCCCCCC)C(=O)Cl.C(CCCCCCCCCCCCCC)C(=O)Cl>>C(COC(=O)CCCCCCCCC)(COC(=O)CCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC"
    # m20 = "[CH2:1]([CH2:3][CH2:6][CH2:7][CH2:8][CH2:9][CH2:10][CH2:11][CH2:12][CH2:13][CH2:14][CH2:15][CH2:16][CH2:17][CH3:18])[C:2](=[O:5])Cl.[CH2:19]([CH2:21][CH2:24][CH2:25][CH2:26][CH2:27][CH2:28][CH2:29][CH2:30][CH2:31][CH2:32][CH2:33][CH2:34][CH2:35][CH3:36])[C:20](=[O:23])Cl>>C(COC(=O)CCC[CH2:31][CH2:32][CH2:33][CH2:34][CH2:35][CH3:36])(CO[C:20](=[O:23])[CH2:19][CH2:21][CH2:24][CH2:25][CH2:26][CH2:27][CH2:28][CH2:29][CH3:30])O[C:2](=[O:5])[CH2:1][CH2:3][CH2:6][CH2:7][CH2:8][CH2:9][CH2:10][CH2:11][CH2:12][CH2:13][CH2:14][CH2:15][CH2:16][CH2:17][CH3:18]"
    # s21 = "C(CCCCCCCC)CCCCCCCBr.C(CCCCCCCC)CCCCCCCBr>>C(CCCCCCCCCCCCCCCC)(CCCCCCCCCCCCCCCC)(C(=O)OCC)C(=O)OCC"
    # m21 = "[CH2:1]([CH2:3][CH2:5][CH2:7][CH2:9][CH2:11][CH2:13][CH2:15][CH3:17])[CH2:2][CH2:4][CH2:6][CH2:8][CH2:10][CH2:12][CH2:14]Br.[CH2:18]([CH2:20][CH2:22][CH2:24][CH2:26][CH2:28][CH2:30][CH2:32][CH3:34])[CH2:19][CH2:21][CH2:23][CH2:25][CH2:27][CH2:29][CH2:31]Br>>C([CH2:31][CH2:29][CH2:27][CH2:25][CH2:23][CH2:21][CH2:19][CH2:18][CH2:20][CH2:22][CH2:24][CH2:26][CH2:28][CH2:30][CH2:32][CH3:34])([CH2:14][CH2:12][CH2:10][CH2:8][CH2:6][CH2:4][CH2:2][CH2:1][CH2:3][CH2:5][CH2:7][CH2:9][CH2:11][CH2:13][CH2:15][CH3:17])(C(=O)OCC)C(=O)OCC"
    # s22 = "C(CCCCCCCC)CCCCCCCI.C(CCCCCCCC)CCCCCCCI>>C(CCCCCCCCCCC)(CCCCCCCCCCCCCCCC)(C(=O)OCC)C(=O)OCC"
    # m22 = "[CH2:1]([CH2:3][CH2:5][CH2:7][CH2:9][CH2:11][CH2:13][CH2:15][CH3:17])[CH2:2][CH2:4][CH2:6][CH2:8][CH2:10][CH2:12][CH2:14]I.[CH2:18]([CH2:20][CH2:22][CH2:24][CH2:26][CH2:28][CH2:30][CH2:32][CH3:34])[CH2:19][CH2:21][CH2:23][CH2:25][CH2:27][CH2:29][CH2:35]I>>[C:23]([CH2:21][CH2:19][CH2:18][CH2:20][CH2:22][CH2:24][CH2:26][CH2:28][CH2:30][CH2:32][CH3:34])([CH2:14][CH2:12][CH2:10][CH2:8][CH2:6][CH2:4][CH2:2][CH2:1][CH2:3][CH2:5][CH2:7][CH2:9][CH2:11][CH2:13][CH2:15][CH3:17])(C(=O)O[CH2:35]C)[C:25](=O)O[CH2:27][CH3:29]"
    # s23 = "C(CCCCCCCCC)CCCCCCCCCl.C(CCCCCCCCC)CCCCCCCCCl>>C(=S)(SCCCCCCCCCCCCCCCCCC)SCCCCCCCCCCCCCCCCCC"
    # m23 = "[CH2:1]([CH2:3][CH2:5][CH2:7][CH2:9][CH2:11][CH2:13][CH2:15][CH2:17][CH3:19])[CH2:2][CH2:4][CH2:6][CH2:8][CH2:10][CH2:12][CH2:14][CH2:16]Cl.[CH2:20]([CH2:22][CH2:24][CH2:26][CH2:28][CH2:30][CH2:32][CH2:34][CH2:36][CH3:38])[CH2:21][CH2:23][CH2:25][CH2:27][CH2:29][CH2:31][CH2:33][CH2:35]Cl>>C(=S)(S[CH2:35][CH2:33][CH2:31][CH2:29][CH2:27][CH2:25][CH2:23][CH2:21][CH2:20][CH2:22][CH2:24][CH2:26][CH2:28][CH2:30][CH2:32][CH2:34][CH2:36][CH3:38])S[CH2:16][CH2:14][CH2:12][CH2:10][CH2:8][CH2:6][CH2:4][CH2:2][CH2:1][CH2:3][CH2:5][CH2:7][CH2:9][CH2:11][CH2:13][CH2:15][CH2:17][CH3:19]"
    # s24 = "C(=O)(OCC)CC(=O)CCC.C(=O)(OCC)CC(=O)CCC>>C(C(=O)CCC)(C(=O)OCC)C(C)C(C(=O)CCC)C(=O)OCC"
    # m24 = "[C:1](=[O:4])([O:3][CH2:6][CH3:9])[CH2:2][C:5](=[O:8])[CH2:7][CH2:10][CH3:11].[C:12](=[O:15])([O:14][CH2:17][CH3:20])[CH2:13][C:16](=[O:19])[CH2:18][CH2:21][CH3:22]>>[CH:2]([C:5](=[O:8])[CH2:7][CH2:10][CH3:11])([C:1](=[O:4])[O:3][CH2:6][CH3:9])C(C)[CH:13]([C:16](=[O:19])[CH2:18][CH2:21][CH3:22])[C:12](=[O:15])[O:14][CH2:17][CH3:20]"
    # s25 = "[C@]12(C)[C@@]([H])(CC[C@@H]1C(C)=O)[C@@]3([H])[C@@]([H])([C@]4(C)C(C[C@@H](O)CC4)=CC3)CC2.C(=O)(O)O.NC>>[C@]12(C)[C@@]([H])(CC[C@]1([H])C(C)=O)[C@@]3([H])[C@@]([H])([C@]4(C)C(=CC(=O)CC4)CC3)CC2"
    # m25 = "[C@:1]12([CH3:5])[C@@:2]([H:24])([CH2:7][CH2:8][C@@H:3]1[C:9]([CH3:14])=[O:13])[C@@:6]3([H:25])[C@@:11]([H:26])([C@:15]4([CH3:19])[C:17]([CH2:20][C@@H:22]([OH:23])[CH2:21][CH2:18]4)=[CH:16][CH2:12]3)[CH2:10][CH2:4]2.C(=O)(O)O.NC>>[C@:1]12([CH3:5])[C@@:2]([H:24])([CH2:7][CH2:8][C@:3]1([H])[C:9]([CH3:14])=[O:13])[C@@:6]3([H:25])[C@@:11]([H:26])([C@:15]4([CH3:19])[C:17](=[CH:20][C:22](=[O:23])[CH2:21][CH2:18]4)[CH2:16][CH2:12]3)[CH2:10][CH2:4]2"
    # s26 = "C(Cl)(Cl)Cl.[C@@]12([H])[C@](C)(CC[C@]([H])([C@]3(C)C(C[C@@H](O)CC3)=CC4)[C@]14[H])[C@@]([H])(C([H])(C)CCCC(C)C)CC2>>[C@@]12([C@]([H])(C[C@@]([H])([C@@]3([H])[C@@](C)([C@@]([H])([C@]([H])(C)CCCC(C)C)CC3)CC4)[C@]45[H])O1)[C@]5(C)CC[C@H](O)C2"
    # m26 = "C(Cl)(Cl)Cl.[C@@:5]12([H:33])[C@:6]([CH3:11])([CH2:10][CH2:16][C@:12]([H:36])([C@:17]3([CH3:23])[C:21]([CH2:25][C@@H:28]([OH:30])[CH2:26][CH2:22]3)=[CH:18][CH2:13]4)[C@:7]14[H:34])[C@@:9]([H:35])([C:15]([H:37])([CH3:20])[CH2:19][CH2:24][CH2:27][CH:29]([CH3:32])[CH3:31])[CH2:14][CH2:8]2>>[C@@:21]12([C@:18]([H])([CH2:13][C@@:7]([H:34])([C@@:5]3([H:33])[C@@:6]([CH3:11])([C@@:9]([H:35])([C@:15]([H:37])([CH3:20])[CH2:19][CH2:24][CH2:27][CH:29]([CH3:32])[CH3:31])[CH2:14][CH2:8]3)[CH2:10][CH2:16]4)[C@:12]45[H:36])O1)[C@:17]5([CH3:23])[CH2:22][CH2:26][C@H:28]([OH:30])[CH2:25]2"
    # result = ""
    # for s, m in (
    #         (s1, m1), (s2, m2), (s3, m3), (s4, m4), (s5, m5), (s6, m6), (s7, m7), (s8, m8), (s9, m9), (s10, m10),
    #         (s11, m11),
    #         (s12, m12), (s13, m13), (s14, m14), (s15, m15), (s16, m16), (s17, m17), (s18, m18), (s19, m19), (s20, m20),
    #         (s21, m21), (s22, m22), (s23, m23), (s24, m24)):
    #     try:
    #         d = get_dict(s, m)
    #         result += json.dumps(d, indent=4, sort_keys=True) + "\n"
    #     except AssertionError:
    #         pass
    #     except TypeError:
    #         pass
    # print(result)

    # smiles = "C(=O)(NCCCCCCCCCC)NCCCCCCCCCC"
    # M, mol = get_adjacent_matrix(smiles)
    # print(connect_shortest_path([12, 13], M))


def main2():
    get_rule("N(CCCC)CCCC.C#C>>N(CCCC)(CCCC)C(=O)C=C",
             "[NH:1]([CH2:3][CH2:5][CH2:7][CH3:9])[CH2:2][CH2:4][CH2:6][CH3:8].[CH:10]#[CH:11]>>[N:1]([CH2:3][CH2:5][CH2:7][CH3:9])([CH2:2][CH2:4][CH2:6][CH3:8])C(=O)[CH:10]=[CH2:11]")


def write_line(db, cursor, data):
    insert_sql1 = "INSERT INTO reactions(r1_fingerprint,r2_fingerprint,p1_fingerprint,r1_smiles,r2_smiles,p1_smiles,d1_reaction,d2_reaction,d3_reaction," \
                  "d1r1_fingerprint,d1r2_fingerprint,d1p1_fingerprint,d1r1_smiles,d1r2_smiles,d1p1_smiles," \
                  "d2r1_fingerprint,d2r2_fingerprint,d2p1_fingerprint,d2r1_smiles,d2r2_smiles,d2p1_smiles," \
                  "d3r1_fingerprint,d3r2_fingerprint,d3p1_fingerprint,d3r1_smiles,d3r2_smiles,d3p1_smiles," \
                  "raw_reaction) " \
                  "VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s);"
    cursor.execute(insert_sql1, data)
    db.commit()
    lastid = cursor.lastrowid
    insert_sql2 = "INSERT INTO decompositions(smiles,fingerprint,reaction_id) VALUES (%s,%s,%s)"
    cursor.execute(insert_sql2, (data[5], data[2], lastid))
    db.commit()
    print("Write Line")


def parse_db_reactions():
    mydb = mysql.connector.connect(
        host="192.168.1.11",  # 数据库主机地址
        user="ChemAdmin",  # 数据库用户名
        passwd="Chem@2018",  # 数据库密码
        database="Chemical"
    )
    mycursor = mydb.cursor()

    mycursor.execute("SELECT * FROM Chemical.reactions")

    results = mycursor.fetchall()  # fetchall() 获取所有记录

    for result in results:
        id = result[0]
        rule1 = result[7]
        rule2 = result[8]
        rule3 = result[9]
        parts = rule.split(">>")
        lhs = parts[0]
        rhs = parts[1]
        reactants = lhs.split(".")
        products = rhs.split(".")
        r1, r2 = reactants[0], reactants[1]
        p1 = products[0]
        r1_fp = smiles2fp(r1)
        assert r1_fp is not None
        r2_fp = smiles2fp(r2)
        assert r2_fp is not None
        p1_fp = smiles2fp(p1)
        assert p1_fp is not None
        if r1_fp > r2_fp:
            r1_fp, r2_fp = r2_fp, r1_fp
            r1, r2 = r2, r1
        if (r1_fp, r2_fp, p1_fp) not in record_map.keys():
            record_map.add((r1_fp, r2_fp, p1_fp))


def write_db():
    def get_fingerprints(rule):
        parts = rule.split(">>")
        lhs = parts[0]
        rhs = parts[1]
        reactants = lhs.split(".")
        products = rhs.split(".")
        r1, r2 = reactants[0], reactants[1]
        p1 = products[0]
        r1_fp = smiles2fp(r1)
        assert r1_fp is not None
        r2_fp = smiles2fp(r2)
        assert r2_fp is not None
        p1_fp = smiles2fp(p1)
        assert p1_fp is not None
        if r1_fp > r2_fp:
            r1_fp, r2_fp = r2_fp, r1_fp
            r1, r2 = r2, r1
        return r1, r2, p1, r1_fp, r2_fp, p1_fp

    mydb = mysql.connector.connect(
        host="192.168.1.11",  # 数据库主机地址
        user="ChemAdmin",  # 数据库用户名
        passwd="Chem@2018",  # 数据库密码
        database="Chemical"
    )
    mycursor = mydb.cursor()
    data = read_all_data()
    count = 1
    num = 1
    NFE, NFC = 0, 0
    AE, TE, VE, OE, SE = 0, 0, 0, 0, 0
    record_set = set()
    for l in data:
        error_msg = ""
        oe = ""
        i = 1
        for item in l:
            try:
                osmi, smi = item[0], item[1]
                if not check_rule(osmi):
                    raise NotFormalException
                r1, r2, p1, r1_fp, r2_fp, p1_fp = get_fingerprints(osmi)
                rule1 = get_rule(osmi, smi, dist=1)[0]
                assert rule1 is not None
                rule2 = get_rule(osmi, smi, dist=2)[0]
                assert rule2 is not None
                rule3 = get_rule(osmi, smi, dist=3)[0]
                assert rule3 is not None
                if (r1_fp, r2_fp) in record_set:
                    raise SkipException
                d1r1, d1r2, d1p1, d1r1_fp, d1r2_fp, d1p1_fp = get_fingerprints(rule1)
                d2r1, d2r2, d2p1, d2r1_fp, d2r2_fp, d2p1_fp = get_fingerprints(rule2)
                d3r1, d3r2, d3p1, d3r1_fp, d3r2_fp, d3p1_fp = get_fingerprints(rule3)
                record_set.add((r1_fp, r2_fp, p1_fp))
                value = (r1_fp, r2_fp, p1_fp, r1, r2, p1, rule1, rule2, rule3,
                         d1r1_fp, d1r2_fp, d1p1_fp,d1r1, d1r2, d1p1,
                         d2r1_fp, d2r2_fp, d2p1_fp,d2r1, d2r2, d2p1,
                         d3r1_fp, d3r2_fp, d3p1_fp, d3r1, d3r2, d3p1,osmi)
                write_line(mydb, mycursor, value)
            except NotFormalException:
                NFE = NFE + 1
            except NotFormalConvertion:
                NFC = NFC + 1
            except AssertionError as e:
                AE = AE + 1
                error_msg = error_msg + (str(item) + "\n" + str(e) + "\n")
            except TypeError as e:
                TE = TE + 1
                error_msg = error_msg + (str(item) + "\n" + str(e) + "\n")
            except ValueError as e:
                VE = VE + 1
                error_msg = error_msg + (str(item) + "\n" + str(e) + "\n")
            except SkipException as e:
                SE = SE + 1
            except Exception as e:
                OE = OE + 1
                oe = oe + (str(item) + "\n" + str(e) + "\n")
            finally:
                print("{:d}_{:d}".format(count, i))
                i = i + 1
                num = num + 1
        count = count + 1
    print("Num:{:d},NotFormalException:{:d},NotFormalConvertion:{:d},AssertionError:{:d},TypeError:{:d},ValueError:{:d},OtherError:{:d},Skip:{:d}".format(num, NFE, NFC, AE, TE, VE, OE, SE))
    mycursor.close()
    mydb.close()


if __name__ == "__main__":
    # write_db()
    # m = Chem.MolFromSmiles('cc1(C#N)c(N)[nH]nc1CC#N')
    # patt = Chem.MolFromSmarts('CC(=NN)C')
    # m = Chem.MolFromSmiles('C(CN)N(C)C')
    # patt = Chem.MolFromSmarts('C(CN)N')
    m = Chem.MolFromSmiles('c(c(ccc1)ccc2O)(c1\C(=N/O)\c(cccc3)c3C4=O)c24')
    patt = Chem.MolFromSmarts('NN')

    print(m.GetSubstructMatches(patt))