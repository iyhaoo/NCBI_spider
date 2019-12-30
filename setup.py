from setuptools import setup, find_packages

requirements = [
    #"numpy>=1.14.0",
    #"pandas>=0.21.0",
    #"tensorflow>=1.13.1,<2.0.0",
    #"h5py>=2.9.0"
]

setup(
    name="ncbispider",
    version="0.0.1",
    author="Hao Yuan",
    author_email="904469382@qq.com",
    description="",
    install_requires=requirements,
    url="https://github.com/iyhaoo/NCBI_spider",
    packages=find_packages(),
    license='Apache License 2.0',
    classifiers=[
        #"Programming Language :: Python :: 3.6",
        #"License :: OSI Approved :: Apache Software License",
        #"Topic :: Scientific/Engineering :: Artificial Intelligence"
    ],
    entry_points={
        'console_scripts': [
            'ncbispider = ncbispider.__main__:main'
        ]},
    #python_requires='>=3.6',
)




