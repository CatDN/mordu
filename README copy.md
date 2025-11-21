# hello_pypi-cat
Simple python package configuration for publishing to pypi (test pypi in this case).

## How to:

see: https://packaging.python.org/en/latest/tutorials/packaging-projects/ for full walkthrough

1. create base package structure (source code in package named directory within the src directory)

2. develop the package within the src directory (include examples)

3. create a tests folder and populate it

4. create metadata files:

    - .gitignore

    - LICENSE

    - pyproject.toml

    - README.md

5. populate meta data

    - anything not meant for public eyes into gitignore

    - choose the license and copy paste to LICENSE (see https://choosealicense.com/)

    - use this example as base for the .toml file

    - add any details to README for the future users

6. generate distribution archives: `py -m build`

7. upload the distribution packages to pypi: `py -m twine upload dist/*
`

    - ensure the current package name is acceptable before attempting