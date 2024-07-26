
# Altbc_analyzer

## 项目简介
ALTBC Analyzer 是一个用于计算硫族相变存储材料晶体结构中角度限制的三体关联函数（the angular-limited three-body correlation, ALTBC）的工具。
该项目使用 ASE 库读取结构文件，并使用近邻列表线性标度算法降低计算量。

## 文件结构
```
project/
├── altbc_analyzer
│   ├── altbc_analyzer.py
│   ├── __init__.py
│   └── neighbor_list.py
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

## 安装步骤
1. 克隆本项目到本地：
   ```bash
   git clone https://github.com/wangchr1617/altbc_analyzer.git
   cd altbc_analyzer
   ```

2. 安装项目依赖：
   ```bash
   pip install -r requirements.txt
   ```

## 使用方法
### 分析和绘图
在项目根目录下运行`analyze_and_plot.py`文件来进行分析和绘图。以下是主要的函数和它们的功能：

- `analyze_file(filename, filetype='POSCAR', cutoff=4.0, theta_min=155, theta_max=180, frame_interval=10)`：分析指定的文件类型，并计算结构文件的ALTBC。
- `plot_results(df, plot_type='scatter')`：根据分析结果绘制散点图或核密度图。

### 示例
以下是一个运行示例：
```bash
python analyze_and_plot.py
```

该命令将分析指定的晶体结构文件并生成相应的绘图结果，保存为`altbc.png`。

## 示例代码片段
```python
from altbc_analyzer.altbc_analyzer import ALTBC_Analyzer
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
欢迎提交 Issue 和 Pull Request 来改进本项目。如果有任何问题或建议，请联系项目维护者[邮箱](wangchr1617@gmail.com)。

## 许可证
本项目基于MIT许可证开源，详情请参阅LICENSE文件。
