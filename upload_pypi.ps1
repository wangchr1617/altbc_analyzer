
Remove-Item -Recurse -Force .\build
Remove-Item -Recurse -Force .\dist
Remove-Item -Recurse -Force .\altbc_analyzer.egg-info

python setup.py sdist bdist_wheel

$env:TWINE_USERNAME = "__token__"
$env:TWINE_PASSWORD = "pypi-AgEIcHlwaS5vcmcCJDgyZDU5YTM1LWIwODEtNDNiNi04MTA1LWQ5NzI3MzE3MGI0NwACKlszLCI4NDcxNThjMS1lODIwLTRjMTUtOTYxYS1lNTk4NGM0NzRkYzEiXQAABiABqDEw8kEEiMj6xguv6ZiJwjJ_A05DhiNL-85rdMdaRQ"

python -m twine upload dist/*
