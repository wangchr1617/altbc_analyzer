
# Altbc_analyzer

## Project Overview
ALTBC Analyzer is a Python library for calculating the angular-limited three-body correlation (ALTBC) function in the crystal structure of chalcogenide phase change materials.
This project uses the ASE library to read structure files and employs a linear scaling algorithm with neighbor lists to reduce computation.

## File Structure
```
project/
├── altbc_analyzer
│   ├── __init__.py
│   ├── altbc_analyzer.py
│   └── neighbor_list.py
├── example
│   └── POSCAR
│   └── XDATCAR
├── LICENSE
├── README.md
├── requirements.txt
├── setup.py
└── Test.ipynb
```

## Installation

```bash
git clone https://github.com/wangchr1617/altbc_analyzer.git
cd altbc_analyzer
pip install -r requirements.txt
```

## Contribution
You are welcome to submit Issues and Pull Requests to improve this project. If you have any questions or suggestions, please contact the project maintainer via [email](mailto:wangchr1617@gmail.com).

## License
This project is open-sourced under the MIT license. For details, please refer to the LICENSE file.
