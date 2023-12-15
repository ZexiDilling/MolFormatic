from rdkit import Chem


file = "C:/Users/phch/Desktop/more_data_files/Prestwick_Chemical_Library_Ver19_384.sdf"
output_file = "C:/Users/phch/Desktop/more_data_files/Prestwick_Chemical_Library_Ver19_384_old_format.sdf"
supplier = Chem.SDMolSupplier(file)
writer = Chem.SDWriter(output_file)
counter = 0
for molecule in supplier:
    counter += 1
    # print(counter)
    try:
        smiles = Chem.MolToSmiles(molecule)
    except:
        print(temp_mol.GetProp("Compound_identifying_number"))
    else:
        # print(counter)
        prop_list = list((molecule.GetPropNames()))
        # print(smiles)
        # print(molecule.GetProp("Compound_identifying_number"))
        # print(counter)
        #
        mol = Chem.MolFromSmiles(smiles)

        for header in prop_list:
            temp_prop = molecule.GetProp(header)

            # print(temp_prop)
            mol.SetProp(header, temp_prop)

        writer.write(mol)
        temp_mol = molecule

writer.close()