.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

The Code
========

Our Git repository contains the Axom source code, documentation, test
suites and all files and scripts used for configuring and building the code.
The repository lives in our
`Github repository <https://github.com/LLNL/axom>`_.

We use `Github <https://github.com/LLNL/axom/issues>`_ for
issue tracking. Please report issues, feature requests, etc. there or send
email to the Axom development team.


Getting the Code
----------------

Access to our repository and Atlassian tools requires membership in the LC
group ``axom``. If you're not in the group, please send email to
'axom-dev@llnl.gov' and request to be added.

SSH keys
^^^^^^^^

If you have not used Github before, you should start by creating and adding your SSH keys to Github. 
Github provides a `good tutorial <https://help.github.com/en/enterprise/2.18/user/github/authenticating-to-github/adding-a-new-ssh-key-to-your-github-account>`_.
Performing these two simple steps will make it easier for you to interact with 
our Git repository without having to repeatedly enter login credentials. 

Cloning the repo
^^^^^^^^^^^^^^^^

To clone the repo into your local working space, type the following::

  $ git clone --recursive git@github.com:LLNL/axom.git

Important notes:

  * You don't need to remember the URL for the Axom repo above. It can be
    found by going to the Axom repo on our Github project and
    clicking on the 'Clone or download' button that is on the upper right of the Axom Github
    page.
  * The ``--recursive`` argument above is needed to pull in Axom's submodules.
    This includes our data directory, which is used for testing, as well as our 
    build system called *BLT*, a standalone product that lives in
    `its own repository <https://github.com/llnl/blt>`_.
    Documentation for BLT can be found `here <https://llnl-blt.readthedocs.io/en/latest/>`_.
  * If you forget to pass the ``--recursive`` argument to the git clone command,
    the following commands can be typed after cloning:

    .. code:: bash

      $ cd axom
      $ git submodule init
      $ git submodule update


Repository Layout
-----------------

If you need to look through the repository, this section explains how it is
organized. If you do not need this information and just want to build the
code, please continue on to the next section.

The top-level Axom directory contains the following directories:

  data
    The optional `axom_data` submodule is cloned here.
    
  host-configs
    Detailed configuration information for platforms and compilers we support.

    See :ref:`hostconfig-label` for more information.
  scripts
    Scripts that we maintain to simplify development and usage tasks
  src
    The bulk of the repo contents.

    Within the **src** directory, you will find the following directories:

    axom
      Directories for individual Axom components (see below)
    cmake
      Axom's build system lives here.

      The BLT submodule is cloned into the **blt** subdirectory.
    docs
      General Axom documentation files
    examples
      Example programs that utilize Axom in their build systems
    thirdparty
      Built-in third party libraries with tests to ensure they are built properly.

In the **axom** directory, you will find a directory for each of the
Axom components. Although there are dependencies among them, each is
developed and maintained in a largely self-contained fashion. Axom
component dependencies are essentially treated as library dependencies.
Each component directory contains subdirectories for the component header
and implementation files, as well as user documentation, examples and tests.
