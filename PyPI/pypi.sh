rm -rf build/ dist/ eeisp.egg-info/
python -m build
#twine upload --repository testpypi dist/*
twine upload --repository pypi dist/*
