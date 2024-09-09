import chemutil

mols = chemutil.MolIterator()
mols.ReadXYZFile('benzene.xyz')
print('File contains {} geometries'.format(mols.Size()))
if mols.Size() > 0:
   mol = mols[0]
   print('Molecule contains {} atoms'.format(mol.NAtoms()))
