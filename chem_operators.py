from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.AtomPairs import Pairs


def get_mol(smiles):
    """
    Gets molecule weight based on smiles code

    :param smiles: Smiles code for compound
    :type smiles: str
    :return: molecule weight in mols
    :rtype: Chem.MolFromSmiles
    """
    return Chem.MolFromSmiles(smiles)


def mw_from_smiles(smiles):
    """
    Gets molecule weight based on smiles code

    :param smiles: Smiles code for compound
    :type smiles: str
    :return: molecule weight in mols
    :rtype: int
    """
    mol = Chem.MolFromSmiles(smiles)
    mw = round(Descriptors.MolWt(mol), 4)
    return mw


def png_string(smiles, size=(250, 100)):
    """
    Translate smiles code to a png-string that can be drawn, to get a 2d drawing of compounds from a smiles code

    :param smiles: Smiles code for a compound
    :type smiles: str
    :param size: The size of the picture
    :type size: (int, int)
    :return: A png in string formate
    :rtype: str
    """
    mol = get_mol(smiles)
    data = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    data.DrawMolecule(mol)
    data.FinishDrawing()
    return data.GetDrawingText()


def _threshold_search(methode, smiles_search, compound_data, threshold, morgan_values):
    sub_search = sub_search_method(methode, morgan_values)
    if not sub_search:
        return None, None
    smiles_1 = smiles_search
    compounds_to_keep = {}
    hit_list = set()

    for compound in compound_data:
        smiles_2 = compound_data[compound]["smiles"]
        smiles_temp = (smiles_1, smiles_2)
        score = sub_search(smiles_temp)

        if type(score) == tuple:
            check_score = score[-1]
        else:
            check_score = score
        if threshold:
            if check_score >= threshold:
                compound_data[compound]["match_score"] = round(check_score, 2)
                hit_list.add(compound)
                compounds_to_keep[compound] = compound_data[compound]

        else:
            compound_data[compound]["match_score"] = round(check_score, 2)
            compounds_to_keep[compound] = compound_data[compound]

    return compound_data, compounds_to_keep


def _clean_up_smiles(smiles_list, unspecified):
    chem_patterns = []
    for smiles in smiles_list:
        temp_chem = Chem.MolFromSmiles(smiles)
        if unspecified:
            # temp_smiles = Chem.MolToSmiles(temp_chem).replace("-", "~")
            temp_smiles = Chem.MolToSmarts(temp_chem)
            temp_chem = Chem.MolFromSmarts(temp_smiles)

        chem_patterns.append(temp_chem)

    return chem_patterns


def _skeleton_search(compound_data, smiles_list, unspecified, testing):
    hit_list = set()
    compounds_to_keep = {}
    chem_patterns = _clean_up_smiles(smiles_list, unspecified)
    for compound_index, compound in enumerate(compound_data):
        if testing:
            temp_chem = Chem.MolFromSmiles(compound)
        else:
            test_compound = compound_data[compound]["smiles"]
            temp_chem = Chem.MolFromSmiles(test_compound)
        for pattern_index, pattern in enumerate(chem_patterns):
            if temp_chem.HasSubstructMatch(pattern):
                hit_list.add(compound)

        if compound in hit_list:
            compounds_to_keep[compound] = compound_data[compound]

    return hit_list, compounds_to_keep


def structure_search(method, threshold, compound_data, smiles_target_list,
                     morgan_values=None, unspecified=True, data_in_list=False):
    """
    Compare molecules with a main smiles code, to see how similar they are.

    :param method: What structure search methode to use. (general, morgan or dice)
    :type method: str
    :param threshold: Minimum similarity score the compound needs
    :type threshold: int
    :param compound_data: Rows from the database with compounds.
    :type compound_data: dict
    :param smiles_target_list: A list of smiles code that are used for targeting the search
    :type smiles_target_list: list
    :param morgan_values: Morgan search values, for determining different variables of the search criteria.
    :type morgan_values: list or None
    :param unspecified: This should ignore angels for bindings for the compound if set to true. (Not sure that it works)        #ToDo Check if this is the case
    :type unspecified: bool
    :param data_in_list: if compound data is in a list instead of a dict formate
    :type data_in_list: bool
    :return: Rows from the database, with compounds under the threshold removed.
    :rtype: dict
    """

    if method.casefold() == "skeleton":
        hit_list, compounds_to_keep = _skeleton_search(compound_data, smiles_target_list, unspecified, data_in_list)
    else:
        hit_list, compounds_to_keep = _threshold_search(method, smiles_target_list[0], compound_data, threshold,
                                                        morgan_values)

    return compounds_to_keep, hit_list


def sub_search_method(methode, morgan_values):
    """
    set the methode to use

    :param methode: What methode to use
    :type methode: str
    :return: Sub_search with the right set-up
    :rtype: function
    """
    if methode.casefold() == "Finger Print":
        sub_search = match_sub_structure_general
    elif methode.casefold() == "morgan":
        sub_search = match_morgan(morgan_values)
    elif methode.casefold() == "dice":
        sub_search = match_dice
    else:
        return None

    return sub_search


def match_sub_structure_general(smiles):
    """
    Uses Fingerprint search to find out how similar two compounds are

    :param smiles: The two smiles codes that are compared
    :type smiles: tuple
    :return: Match score for how similar the two compounds are
    :rtype: int
    """
    mols = [get_mol(smile) for smile in smiles]
    fp1, fp2 = [FingerprintMols.FingerprintMol(mol) for mol in mols]

    match_score = DataStructs.FingerprintSimilarity(fp1, fp2) * 100
    return match_score


def match_morgan(smiles, bonds=None, n_bits=None, chirality=None, features=None):
    """
    Using Morgan search to find out how similar two compounds are
    NEEDS TO READ UP ON THIS! ! !

    :param smiles: The two smiles codes that are compared
    :type smiles: tuple
    :param bonds:
    :type bonds:
    :param n_bits:
    :type n_bits: int
    :param chirality:
    :type chirality: bool
    :param features:
    :type features: bool
    :return: Match score for how similar the two compounds are
    :rtype: int
    """
    mols = [get_mol(smile) for smile in smiles]
    fp1, fp2 = [AllChem.GetMorganFingerprintAsBitVect(mol, 2)
                for mol in mols]
    # fp1, fp2 = [AllChem.GetMorganFingerprintAsBitVect(mol, radius=bonds, nBits=n_bits,
    #                                                   useChirality=chirality, useFeatures=features)
    #             for mol in mols]

    match_score_1 = DataStructs.FingerprintSimilarity(fp1, fp2) * 100
    match_score_2 = DataStructs.TanimotoSimilarity(fp1, fp2) * 100
    return match_score_1, match_score_2


def match_dice(smiles):
    """
    uses Dice search to find out how similar two compounds are

    :param smiles: The two smiles codes that are compared
    :type smiles: tuple
    :return: Match score for how similar the two compounds are
    :rtype: int
    """
    mols = [get_mol(smile) for smile in smiles]
    ap1, ap2 = [Pairs.GetAtomPairFingerprint(mol) for mol in mols]
    dice_match_score = DataStructs.DiceSimilarity(ap1, ap2) * 100
    return dice_match_score


if __name__ =="__main__":

    # #smiles1 = "c1nccc2n1ccc2"
    # smiles1 = "CCO"
    # smiles2 = "CCOC"
    # smiles = (smiles2, smiles1)
    #
    # chirality = True
    # features = False
    # n_bits = 2048
    # bound_range = 2
    #
    # size = (100, 100)

    smiles_list = ['CC1=NN(CC1)C2=CC=CC=C2', "NC1=CC=C(N=C(S)S2)C2=C1"]

    compound_list = [
        "N#CCNC(=O)C1=NN(c2ccccc2)C(c2ccccc2)C1",
        "Cc1c(CNC(=O)C2CC(C(=O)O)=NN2c2ccccc2)cnn1C",
        "Cc1cc(CNC(=O)C2=NN(c3ccccc3)C(C(N)=O)C2)no1",
        "NC(=O)C1CC(C(=O)Nc2cc[nH]n2)=NN1c1ccccc1",
        "COCC(C)NC(=O)C1=NN(c2ccc(F)cc2)C(C(N)=O)C1",
        "CC1=NN(c2ccc(-c3nnc(-c4cnccn4)o3)cc2)CC1",
        "CC1=NN(c2ccc(C(=O)Nc3ccc(CS(C)(=O)=O)cc3)cc2)CC1",
        "O=C(O)C1=NN(c2ccccc2)C(=O)C1CCO",
        "CC1=NN(c2ccc(C(=O)NCCc3nc4c(s3)CCCC4)cc2)CC1",
        "COc1ccc(NC(C)=O)c(NC(=O)c2ccc(N3CCC(C)=N3)cc2)c1",
        "CC1=NN(c2ccc(C(=O)N(CCC#N)Cc3ccco3)cc2)CC1",
        "CC1=NN(c2ccc(C(=O)N3CCCC(NC(=O)C(C)C)C3)cc2)CC1",
        "O=C(O)C1=NN(c2ccccc2)C(C(=O)NC2CCCOC2)C1",
        "Cc1nc(CNC(=O)c2ccc(N3N=C(Cc4ccccc4)CC3=O)cc2)no1",
        "Nc1ccc2nc(SCC(=O)Nc3ccccc3)sc2c1",
        "Nc1ccc2nc(SCC(=O)NCc3ccco3)sc2c1",
        "Nc1ccc2nc(SCC(=O)Nc3ccccc3F)sc2c1",
]

    compound_data, hit_list, compounds_to_keep = structure_search("skeleton_search", None, compound_list, smiles_list, unspecified=False, data_in_list=True)
    print(len(pattern_list))
    compound_data, hit_list, compounds_to_keep = structure_search("skeleton_search", None, compound_list, smiles_list, unspecified=True, data_in_list=True)
    print(len(pattern_list))

    # pattern = Chem.MolFromSmarts('[#6]12~[#6][#6]~[#6][#6]~[#6]1[#6]3~[#6][#6]~[#6][#6]4~[#6]3[#6]2~[#6][#6]~[#6]4')
