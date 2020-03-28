import setuptools

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="fedorov",  # Replace with your own username
    version="0.0.0",
    author="Pengji Zhou",
    author_email="zhoupj@umich.edu",
    description="A python code repo to initialize different crystal structures "
                "based on Pearson symbol or space group",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/zhou-pj/crystal_database",
    packages=['fedorov'],
    package_data={'fedorov': ['crystal_data/*.csv', 'crystal_data/*.json']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Topic :: Scientific/Engineering",
    ],
    install_requires=[
        'numpy>=1.10',
        'spglib>=1.9.0',
        'pandas>=0.20.0',
        'rowan>=1.0.0'],
    python_requires='>=3.3',
)
