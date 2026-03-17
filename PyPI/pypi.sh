rm -rf build/ dist/ eeisp.egg-info/
#python setup.py bdist_wheel
python -m build
#twine upload --repository testpypi dist/*
twine upload --repository pypi dist/*
