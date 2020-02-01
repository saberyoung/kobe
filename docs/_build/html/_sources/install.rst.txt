Installation
===================================

Requirements
------------

PSTools depends on several python libraries, e.g. `healpy <https://healpy.readthedocs.io/en/latest/>`_, 
`matplotlib <https://matplotlib.org/>`_, `astropy <https://www.astropy.org/>`_, etc
(see `full list <https://github.com/saberyoung/pstools/blob/master/requirements.txt>`_).
If you adopt pip approach, you need to do: pip install pstools --user, 
all the libraries that needed by PSTools can be installed via pip.
If you use source setup installation, you can do: pip install -r requirements.txt,
and pip would install them automatically so you do not need to 
install them yourself.
If you prefer bash to the setup, `astrocond <https://astroconda.readthedocs.io/en/latest/>`_ 
is suggested, and then you need to install a few libraries if needed.

As an option, PSTools employ `pygcn <https://github.com/lpsinger/pygcn>`_ to 
listens for VOEvent XML format notices over
the custom TCP/IP VOEvent Transport Protocol,
Thus, if you intend to activate such function,
please verify your machine has network, which is not private.

Note that only Linux and MAC OS have been tested, for Windows everything reamins unknown.
Considering the Python version, the tool is working for Python 2/3.

Bash installation (RECOMMENDED)
---------------------------------------------

You can download the source files from `pstools repository <https://github.com/saberyoung/pstools.git>`_ with 
the ``git`` command::

    git clone https://github.com/saberyoung/pstools.git
    
which would create one directory appear in the current path, e.g. ``/home/xxx/pstools/``.
Then, you need to add one line in your bash (if you use bash)::

    alias pstools_init="source /home/xxx/pstools/pst.bash"

After you source your bash file, you can then initial pstools by::

    source activate iraf27/astroconda
    pstools_init

This approach would not install the pipeline, but instead call python via envirnmental defination,
which is light and easy.

Source installation with Pypi
-------------------------------

It is possible to build the latest ``pstools`` with `pip <http://www.pip-installer.org>`_ ::

    pip install --user pstools

If you have installed with ``pip``, you can keep your installation up to date
by upgrading from time to time::

    pip install --user --upgrade pstools

Almost-as-quick installation from official source release
----------------------------------------------------------

PSTools is also available in the
`Python Package Index (PyPI) <https://pypi.org/project/pstools/>`_. You can
download it with::

    curl -O https://files.pythonhosted.org/packages/1c/d5/42cb34cd80b1049b2f4352f17e9277344a25ee07edc8557d9abf9e963147/pstools-0.0.3.tar.gz

and build it with::

    tar -xzf pstools-0.0.3.tar.gz    
    python setup.py install --user

Check
-----

If everything goes fine, you can test it::

    python

>>> import pst
>>> pst.__version__

Clean
-----

When you run "python setup.py", temporary build products are placed in the
"build" directory. If you want to clean out and remove the ``build`` directory,
then run::

    python setup.py clean --all

Uninstall
-----------

For uninstallation, one can easily delete the files directly.
In order to know the file path, you should start python in correct environment and do::

>>> import pst
>>> pst.__path__

Another approach is to remove via pip::

    pip uninstall pstools
