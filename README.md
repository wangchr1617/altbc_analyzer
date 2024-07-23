
# altbc_analyzer

## 项目简介
ALTBC Analyzer 是一个用于计算硫族相变存储材料晶体结构中角度限制的三体关联函数（the angular-limited three-body correlation, ALTBC）的工具。
该项目使用 ASE 库读取结构文件，并使用近邻列表线性标度算法降低计算量。

## 文件结构
```
project/
├── altbc_analyzer
│   ├── altbc_analyzer.py
│   ├── __init__.py
│   ├── neighbor_list.py
│   └── __pycache__
│       ├── altbc_analyzer.cpython-39.pyc
│       ├── __init__.cpython-39.pyc
│       └── neighbor_list.cpython-39.pyc
├── analyze_and_plot.py
├── example
│   └── POSCAR
├── LICENSE
├── README.md
├── requirements.txt
├── setup.py
├── upload_pypi.sh
└── upload_pypi_test.sh
```

### 文件说明
- `altbc_analyzer/` 目录下包含了主要的计算和分析逻辑：
  - `altbc_analyzer.py`：定义了`ALTBC_Analyzer`类，用于计算ALTBC。
  - `neighbor_list.py`：定义了`NeighborList`类，用于生成近邻列表。
  - `__init__.py`：使得该目录可以作为Python模块导入。
- `analyze_and_plot.py`：包含文件分析、数据处理和绘图的主函数。
- `requirements.txt`：列出了项目所需的Python库。

## 安装步骤
1. 克隆本项目到本地：
   ```bash
   git clone https://github.com/wangchr1617/altbc_analyzer.git
   cd altbc_analyzer
   ```
2. 创建并激活虚拟环境（可选）：
   ```bash
   python -m venv venv
   source venv/bin/activate  # 在Windows上使用 venv\Scripts\activate
   ```
3. 安装项目依赖：
   ```bash
   pip install -r requirements.txt
   ```

## 使用方法
### 分析和绘图
在项目根目录下运行`analyze_and_plot.py`文件来进行分析和绘图。以下是主要的函数和它们的功能：

- `analyze_file(filename, filetype='POSCAR', cutoff=4.0, theta_min=155, theta_max=180, frame_interval=10)`：分析指定的文件类型，并计算原子三元组的几何性质。
- `plot_results(df, plot_type='scatter')`：根据分析结果绘制散点图或密度图。

### 示例
以下是一个运行示例：
```bash
python analyze_and_plot.py
```

该命令将分析指定的晶体结构文件并生成相应的绘图结果，保存为`altbc.png`。

## 示例代码片段
```python
from altbc_analyzer import ALTBC_Analyzer
from ase.io import read
import matplotlib.pyplot as plt

def analyze_and_plot(filename):
    analyzer = ALTBC_Analyzer(cutoff=4.0, theta_min=155, theta_max=180)
    atoms = read(filename)
    df = analyzer.compute_altbc(atoms)
    if not df.empty:
        plot_results(df, plot_type='scatter')
        plt.savefig('altbc.png', bbox_inches='tight')

if __name__ == "__main__":
    analyze_and_plot('your_structure_file.vasp')
```

## 注意事项
- 如果要分析的结构文件包含超过1000个原子，建议使用POSCAR文件类型。
- 确保输入文件的路径和文件名正确。

## 贡献
欢迎提交Issue和Pull Request来改进本项目。如果有任何问题或建议，请联系项目维护者。

## 许可证
本项目基于MIT许可证开源，详情请参阅LICENSE文件。
