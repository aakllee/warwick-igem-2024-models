import numpy as np 
import matplotlib.pyplot as plt 

class cell:
    def __init__(self, x, y, cell_density=0, chem_density=0, ln_density=0):
        self.x = x
        self.y = y
        self.cell_density = cell_density
        self.chem_density = chem_density 
        self.ln_density = ln_density
    
    def update_cell(self, new_density):
        self.cell_density = new_density
    
    def update_chem(self, new_density):
        self.chem_density = new_density

    def unpdate_ln(self, new_density):
        self.ln_density = new_density

    def info(self):
        print(f"Cell Density: {self.cell_density}, Chemotractor Concentration: {self.chem_density}, Lenthanides Concentration: {self.ln_density}")