from setuptools import setup

with open("requirements.txt") as f:
    install_requires = [line.strip() for line in f if line.strip() and not line.startswith("#")]

setup(
    name="tb-dashboard",
    version="1.0.0",
    author="Fernando Falat",
    author_email="fernandofalat@proton.me",
    description="A web-based genomic explorer for Mycobacterium tuberculosis drug resistance mutations using the WHO catalogue.",
    py_modules=["app", "data_utils", "coordinate_calculator"],
    install_requires=install_requires,
)
