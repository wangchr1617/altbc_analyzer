
Remove-Item -Recurse -Force .\build
Remove-Item -Recurse -Force .\dist
Remove-Item -Recurse -Force .\altbc_analyzer.egg-info

python setup.py sdist bdist_wheel

$env:TWINE_USERNAME = "__token__"
$env:TWINE_PASSWORD = "pypi-AgENdGVzdC5weXBpLm9yZwIkMDk1YjNhYTItZjRkYS00OTdmLTg4M2ItZmNmYmI4MjkzZDRhAAIqWzMsIjFlZmFjNWJlLThhMGQtNGZhYi04NDkwLWFkYzVkYTExYTljOCJdAAAGIMzTr3zRIuhd3PvMhSEJ4Y9oBGkNN2lpO2vXbmMj8V9Q"

python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
