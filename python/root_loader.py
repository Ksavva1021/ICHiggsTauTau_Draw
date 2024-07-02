# root_loader.py

import ROOT
import os

# Path to the shared library
lib_path = os.path.join(os.getcwd(), 'lib', 'libMultiDraw.so')

# Load the shared library
ROOT.gSystem.Load(lib_path)

# Declare the MultiDraw function
ROOT.gInterpreter.Declare('''
        extern void MultiDraw(TTree *inTree,
            TObjArray *Formulae, TObjArray *Weights, TObjArray *Hists,
            UInt_t ListLen);
''')
