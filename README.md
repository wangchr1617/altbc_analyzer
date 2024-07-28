
# Altbc_analyzer

## Project Overview
ALTBC Analyzer is a Python library for calculating the angular-limited three-body correlation (ALTBC) function in the crystal structure of chalcogenide phase change materials.
This project uses the ASE library to read structure files and employs a linear scaling algorithm with neighbor lists to reduce computation.

## File Structure
```
project/
├── altbc_analyzer
│   ├── altbc_analyzer.py
│   ├── __init__.py
│   └── neighbor_list.py
├── analyze_and_plot.py
├── example
│   └── POSCAR
├── LICENSE
├── README.md
├── requirements.txt
└── setup.py
```

## Installation

### By pip 

```shell
$ pip install altbc_analyzer
```

### From Source

```bash
git clone https://github.com/wangchr1617/altbc_analyzer.git
cd altbc_analyzer
pip install -r requirements.txt
```

## Usage
### Analysis and Plotting
Run the `analyze_and_plot.py` file in the project root directory for analysis and plotting. Below are the main functions and their features:

- `analyze_file(filename, filetype='POSCAR', cutoff=4.0, theta_min=155, theta_max=180, frame_interval=10)`: Analyzes the specified file type and calculates the ALTBC of the structure file.
- `plot_results(df, plot_type='scatter')`: Plots scatter or kernel density graphs based on the analysis results.

### Example
Here is an example of running the script:
```bash
python analyze_and_plot.py
```

This command will analyze the specified crystal structure file and generate the corresponding plot, saved as `altbc.png`.

## Example Code Snippet
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

## Notes
- If the structure file to be analyzed contains more than 1000 atoms, it is recommended to use the POSCAR file type.
- Ensure the input file path and name are correct.

## Contribution
You are welcome to submit Issues and Pull Requests to improve this project. If you have any questions or suggestions, please contact the project maintainer via [email](mailto:wangchr1617@gmail.com).

## License
This project is open-sourced under the MIT license. For details, please refer to the LICENSE file.
