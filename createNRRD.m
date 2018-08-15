load('Vb');
filename = 'Data/binModel.nrrd';
pixelspacing = [0.1953, 0.1953, 0.3400];
origin = [0, 0, 0];
encoding = 'raw';

nrrdWriter(filename, Vb, pixelspacing, origin, encoding);