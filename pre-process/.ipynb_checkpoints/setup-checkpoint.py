# -*- coding: utf-8 -*-
def main():
    setup(
        name="SCREE",
        version="0.0.3",
        package_dir = {'SCREE':'SCREE'},
        packages =['SCREE'],
        install_requires=[
            'ruamel.yaml'
        ],
        include_package_data = True,
        scripts = ['SCREE/SCREE.py'],
    )


if __name__ == "__main__":
    try:
        from setuptools import setup, find_packages
        main()
    except ImportError:
        print("Can not load setuptools!")