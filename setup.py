from setuptools import setup, find_packages

setup(
    name="mon_package",
    version="1.1.2",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
)