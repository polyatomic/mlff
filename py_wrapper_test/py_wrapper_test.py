import chemutil

mols = chemutil.MolIterator()
mols.ReadXYZFile('benzene.xyz')
print('File contains {} geometries'.format(mols.Size()))
