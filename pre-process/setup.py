# -*- coding: utf-8 -*-
def main():
    setup(
        name="SCREE",
        version="0.0.4",
        package_dir = {'SCREE':'SCREE'},
        package_data = {'SCREE':['SCREE/*']},
        packages =['SCREE'],
        include_package_data = True,
        scripts = ['SCREE/SCREE.py', 'SCREE/bin.R', 'SCREE/analysis.R'],
        description = "SCREE(Single-cell CRISPR scREen data analyses and pErturbation modeling) is a workflow for single-cell CRISPR screens data processing and analysis.", 
        entry_points={
        'console_scripts': [
            'SCREE = SCREE:main'
        ]
        }
    )


if __name__ == "__main__":
    try:
        from setuptools import setup, find_packages
        main()
    except ImportError:
        print("Can not load setuptools!")