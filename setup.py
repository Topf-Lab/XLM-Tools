import setuptools
import subprocess as sp

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="XLM-Tools", # Replace with your own username
    version="1.0",
    author="Matthew James Sinnott, Joshua Matthew Allen Bullock, Konstantinos Thalassinos, Maya Topf",
    author_email="msinno01@mail.bbk.ac.uk",
    description="Crosslink Modelling Tools (XLM-Tools): a combined tool to score model protein structures using crosslinks an monolinks.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Topf-Lab/XLM-Tools",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: GNU :: GPLv3.0",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)

sp.call(["chmod","+x","alias.bash"])
sp.call(["bash", "./alias.bash"])