NCBI_spider
===========

|PyPI|

.. |PyPI| image:: https://img.shields.io/pypi/v/NCBIspider.svg
    :target: https://pypi.org/project/ncbispider

* Deprecated : the NCBI server seems not support resume so that this repository is useless

Installation
------------

**Install NCBI_spider with pip**
  To install with ``pip``, run the following from a terminal::

    pip install ncbispider

**Install NCBI_spider from GitHub**
  To clone the repository and install manually, run the following from a terminal::

    git clone git://github.com/iyhaoo/NCBI_spider.git

    cd ncbispider

    python setup.py install

Usage
-----
Run NCBI_spider::

  ncbispider --dataset="GSE109762,GSE124557" --out-dir=E:/monkey_single_cell/ChIP
