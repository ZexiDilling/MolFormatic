
from rdkit import Chem


file = "C:/Users/phch/Desktop/more_data_files/Prestwick_Chemical_Library_Ver19_384.sdf"
output_file = "C:/Users/phch/Desktop/more_data_files/Prestwick_Chemical_Library_Ver19_384_old_format.sdf"
supplier = Chem.SDMolSupplier(file)
writer = Chem.SDWriter(output_file)

for molecule in supplier:
    smiles = Chem.MolToSmiles(molecule)
    prop_list = list((molecule.GetPropNames()))

    mol = Chem.MolFromSmiles(smiles)

    for header in prop_list:
        temp_prop = molecule.GetProp(header)
        # print(temp_prop)
        mol.SetProp(header, temp_prop)

    writer.write(mol)

writer.close()


