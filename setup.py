import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="basicsynbio",
    version="0.1.0",
    url="https://github.com/hainesm6/basicsynbio",
    license='MIT',
    author="Matthew Haines",
    author_email="hainesm6@gmail.com",
    description="An open-source Python API to facilitate BASIC DNA Assembly workflows",    
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(exclude=('tests',)),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    install_requires=[
        "biopython>=1.76"
    ],
)