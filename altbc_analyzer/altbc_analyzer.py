from ase.io import read
from scipy.stats import gaussian_kde
from .neighbor_list import NeighborList
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

class ALTBC_Analyzer:
    
    def __init__(self, filename, filetype='SingleSystem', frame_interval=10, center_element=None, xmin=2.5, xmax=3.8, grid_size=0.001, cutoff=4.0, theta_min=155, theta_max=180, apply_mic=False, apply_half=False, savegrid=False):
        """
        初始化 ALTBC_Analyzer 类的实例。
        
        参数：
        filename (str): 输入结构文件，需要包含晶格信息
        filetype (str): 文件类型，仅支持 'SingleSystem' 或 'MultiSystem'
        frame_interval (int): 处理的帧间隔（仅适用于 'MultiSystem' 类型）
        center_element (str): 是否要指定中心原子元素类型
        xmin (float): 网格的最小值
        xmax (float): 网格的最大值
        grid_size (float): 网格大小
        cutoff (float): 原子键合截断距离
        theta_min (float): ALTBC 最小角度（用于过滤）
        theta_max (float): ALTBC 最大角度（用于过滤）
        apply_mic (bool): 是否应用最小图像约束
        apply_half (bool): 是否应用半空间约束
        savegrid (bool): 是否保存画图原始数据
        """
        self.filename = filename
        self.filetype = filetype
        self.frame_interval = frame_interval
        self.center_element = center_element
        self.nframe = 1
        self.xmin = xmin
        self.xmax = xmax
        self.grid_size = grid_size
        self.xbins = int((self.xmax - self.xmin) / self.grid_size)
        self.grid = np.zeros((self.xbins, self.xbins))
        self.cutoff = cutoff
        self.theta_min = theta_min
        self.theta_max = theta_max
        self.apply_mic = apply_mic
        self.apply_half = apply_half
        self.savegrid = savegrid

    def read_structures(self):
        """
        读取文件中的结构数据，支持 'SingleSystem' 和 'MultiSystem' 。
        """
        if self.filetype == 'SingleSystem':
            return [read(self.filename)]
        elif self.filetype == 'MultiSystem':
            frames = read(self.filename, index=slice(0, None, self.frame_interval))
            self.nframe = len(frames)
            return frames
        else:
            print("Error: Invalid filetype. Please specify 'SingleSystem' or 'MultiSystem'.")
            sys.exit(1)

    def generate_neighbor_list(self, atoms):
        """
        生成邻居列表和原子间距。
        """
        nl = NeighborList(self.cutoff, self.apply_mic, self.apply_half)
        return nl.find_neighbors(atoms)

    def find_triplets_with_distances(self, center_atom_index, neighbor_indices, neighbor_distances, i_list, j_list):
        """
        查找所有的三元组及其对应的原子间距。
        """
        triplets = []
        for first_neighbor in neighbor_indices:
            for second_neighbor in neighbor_indices:
                if first_neighbor != second_neighbor:
                    distance_1 = neighbor_distances[(i_list == center_atom_index) & (j_list == first_neighbor)]
                    distance_2 = neighbor_distances[(i_list == center_atom_index) & (j_list == second_neighbor)]
                    triplets.append([first_neighbor, center_atom_index, second_neighbor, distance_1[0], distance_2[0]])
        return np.array(triplets)

    def _compute_altbc(self, atoms):
        """
        计算单个结构的 ALTBC 值。
        """
        i_list, j_list, neighbor_distances = self.generate_neighbor_list(atoms) 
        triplets = []

        if self.center_element is not None:
            center_atom_indices = [i for i, atom in enumerate(atoms) if atom.symbol == self.center_element]
        else:
            center_atom_indices = range(len(atoms))
            
        for center_atom_index in center_atom_indices:
            neighbor_indices = j_list[i_list == center_atom_index]
            if len(neighbor_indices) > 1:
                triplet = self.find_triplets_with_distances(center_atom_index, neighbor_indices, neighbor_distances, i_list, j_list)
                triplets.extend(triplet)

        if not triplets:
            return pd.DataFrame()

        triplets = np.array(triplets)
        angles = atoms.get_angles(triplets[:, :3].astype(int), mic=True)
        
        df = pd.DataFrame(np.hstack([triplets, angles[:, None]]), columns=["A", "B", "C", "AB", "BC", "angle_ABC"])
        df = df[(df["angle_ABC"] >= self.theta_min) & (df["angle_ABC"] <= self.theta_max)].reset_index(drop=True)
        df = df[(df["AB"] <= self.cutoff) & (df["BC"] <= self.cutoff)].reset_index(drop=True)
        return df

    def compute_altbc(self):
        """
        计算所有结构的 ALTBC 值。
        """
        structures = self.read_structures()
        if self.filetype == 'SingleSystem':
            return self._compute_altbc(structures[0])
        elif self.filetype == 'MultiSystem':
            data_frames = [self._compute_altbc(structure) for structure in structures]
            return pd.concat(data_frames, ignore_index=True) if data_frames else pd.DataFrame()

    def _compute_grid(self, df):
        """
        根据给定的 DataFrame 计算并更新网格数据。
        """
        x = np.array(df["AB"])
        y = np.array(df["BC"])
        _grid = np.zeros((self.xbins, self.xbins))
        for i in range(len(x)):
            x_index = int((x[i] - self.xmin) / self.grid_size)
            y_index = int((y[i] - self.xmin) / self.grid_size)
            
            if 0 <= x_index < self.xbins and 0 <= y_index < self.xbins:
                _grid[y_index, x_index] = 1  
        self.grid = self.grid + _grid
        
    def compute_grid(self):
        """
        计算网格，并将结果更新到 `self.grid`。
        """        
        structures = self.read_structures()
        for structure in structures:
            df = self._compute_altbc(structure)
            self._compute_grid(df)
        if self.savegrid:
            np.savetxt("grid.txt", self.grid, fmt = '%f', delimiter = ',')

    def _plot_common_settings(self):
        """
        设置常见的绘图参数，如标签、范围和网格。
        """
        plt.xlabel("Distance r1 (Å)")
        plt.ylabel("Distance r2 (Å)")
        plt.xlim(self.xmin, self.xmax)
        plt.ylim(self.xmin, self.xmax)
        plt.xticks(np.arange(self.xmin+0.1, self.xmax+0.1, 0.2))
        plt.yticks(np.arange(self.xmin+0.1, self.xmax+0.1, 0.2))
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.tight_layout()
    
    def plot_grid(self):
        """
        绘制网格图，并保存为图片。
        """
        self.compute_grid()
        self.grid = self.grid / np.max(self.grid)
        plt.imshow(self.grid, origin='lower', extent=(self.xmin, self.xmax, self.xmin, self.xmax), cmap='jet', interpolation='nearest')
        plt.colorbar(label='P(r1,r2)')
        self._plot_common_settings()

    def plot_scatter(self):
        """
        绘制散点图，并保存为图片
        """
        df = self.compute_altbc()
        x, y = np.array(df["AB"]), np.array(df["BC"])
        scatter = plt.scatter(x, y, s=1)
        self._plot_common_settings()
        