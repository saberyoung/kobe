Welcome to the pstools documentation
=====================================

`Pointing Scheduler Tools (PSTools) <https://sngyang.com/pstools>`_ is a Python package to help develop observing strategies for astronomical telescopes.

`pstools` provides utilities to:

* handle and evaluate trigger (e.g. GW, GRB, or AMON) map via `Healpy/HEALPix <https://healpy.readthedocs.io/en/latest/index.html>`_.
* generate pointing list for specific telescope, depending on its searching strategy, i.e. tiling or galaxy monitoring search.
* optimize observing strategy in order to effectively cover trigger volume, with inclusion of visibility, and the detection efficiency.
* visualize the above process, circulate to users via various of approaches.
  

Installation
----------------

.. toctree::
   :maxdepth: 1

   install

Tutorial
------------

.. toctree::
   :maxdepth: 1

   tutorial
   
Reference
-------------

.. toctree::
   :maxdepth: 2

   pipeline
   view
   cookbook

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

FAQ
------------

.. toctree::
   :maxdepth: 1

   faq

License
---------

.. toctree::
   :maxdepth: 1

   license

Author
-----------
Sheng Yang
http://www.sngyang.com
