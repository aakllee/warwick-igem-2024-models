import numpy as np 
import matplotlib.pyplot as plt
from cell import cell
import os
import imageio

class grid:
    def __init__(self, n, dx, dy, chi, D_chem, D_cell, dt):
        self.dim = n
        self.dx = dx
        self.dy = dy
        self.laplacian_cell = np.zeros((self.dim, self.dim))
        self.laplacian_chem = np.zeros((self.dim, self.dim))
        self.chem_matrix = np.zeros((self.dim, self.dim))
        self.cell_matrix = np.zeros((self.dim, self.dim))
        self.ln_matrix = np.zeros((self.dim, self.dim))
        self.grad_cell = np.zeros((2, self.dim, self.dim))
        self.grad_chem = np.zeros((2, self.dim, self.dim))
        self.D_chem = D_chem
        self.D_cell = D_cell
        self.chi = chi
        self.dt = dt 

        self.cells = []
        for i in range(n):
            celli = []
            for j in range(n):
                celli.append(cell(i,j))
            self.cells.append(celli)
    
    def info(self):
        print(self.cells)
    
    def set_cell(self, cell_matrix):
        self.cell_matrix = np.copy(cell_matrix)
        for i in range(self.dim):
            for j in range(self.dim):
                self.cells[i][j].update_cell(cell_matrix[i][j])

    def set_chem(self, chem_matrix):
        self.chem_matrix = np.copy(chem_matrix)
        for i in range(self.dim):
            for j in range(self.dim):
                self.cells[i][j].update_chem(chem_matrix[i][j])
    
    def set_ln(self, ln_matrix):
        self.ln_matrix = np.copy(ln_matrix)
        for i in range(self.dim):
            for j in range(self.dim):
                self.cells[i][j].update_ln(ln_matrix[i][j])
    
    def calc_laplacian(self):
        self.laplacian_cell = np.copy((np.roll(self.cell_matrix, 1, axis=0) + np.roll(self.cell_matrix, -1, axis=0) + np.roll(self.cell_matrix, 1, axis=1) + np.roll(self.cell_matrix, -1, axis=1) - 4 * self.cell_matrix) / (self.dx * self.dy))
        self.laplacian_chem = np.copy((np.roll(self.chem_matrix, 1, axis=0) + np.roll(self.chem_matrix, -1, axis=0) + np.roll(self.chem_matrix, 1, axis=1) + np.roll(self.chem_matrix, -1, axis=1) - 4 * self.chem_matrix) / (self.dx * self.dy))
    
    def calc_grad(self):
        self.grad_cell[0] = np.copy((np.roll(self.cell_matrix, 1, axis=0) - np.roll(self.cell_matrix, -1, axis=0)) / (2 * self.dx))
        self.grad_cell[1] = np.copy((np.roll(self.cell_matrix, 1, axis=1) - np.roll(self.cell_matrix, -1, axis=1)) / (2 * self.dy))

        self.grad_chem[0] = np.copy((np.roll(self.chem_matrix, 1, axis=0) - np.roll(self.chem_matrix, -1, axis=0)) / (2 * self.dx))
        self.grad_chem[1] = np.copy((np.roll(self.chem_matrix, 1, axis=1) - np.roll(self.chem_matrix, -1, axis=1)) / (2 * self.dy))

    def update(self):
        self.calc_grad()
        self.calc_laplacian()
        new_cell = self.cell_matrix + self.dt * ( self.chi * ( self.grad_cell[0] * self.grad_chem[0] + self.grad_cell[1] * self.grad_chem[1]) + self.laplacian_cell * self.D_cell) 
        new_chem = self.chem_matrix + self.dt * ( self.D_chem * self.laplacian_chem + self.cell_matrix - self.chem_matrix)
        self.set_cell(new_cell)
        self.set_chem(new_chem)

    
    def display(self, x, y):
        self.cells[x][y].info()

newgrid = grid(1000, 1, 1, 1, 1, 1, 0.1)
newgrid.set_cell(np.random.rand(1000, 1000))
newgrid.set_chem(np.random.rand(1000, 1000))
print(newgrid.cell_matrix)

output_dir = "Frames"  # Directory to save frames

# Create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for t in range(1000):
    
    # Save frame
    plt.imshow(newgrid.cell_matrix, cmap='viridis', origin='lower', extent=[0, newgrid.dim * newgrid.dx, 0, newgrid.dim * newgrid.dx])
    plt.colorbar(label='Cell Density')
    plt.title(f"Time step: {t}")
    plt.savefig(f"{output_dir}/frame_{t:04d}.png")
    plt.close()

    newgrid.update()

