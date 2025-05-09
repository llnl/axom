.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _misctasks-label:

********************************
Miscellaneous Development Items
********************************

This section describes various development tasks that need to be 
performed at times and which are not covered in other sections.

===========================
Updating Copyright Headers
===========================

Many files in Axom contain a copyright header that includes dates. These headers
need to be updated in new calendar years. Axom's ``scripts`` directory contains
tools that can update the date across Axom's sources. Before running the tools,
they need to be updated so they contain the new year. This can be done by replacing
the previous year ``(Y-1)`` to the new current year ``(Y)`` and then repeating that
step for the last previous year ``(Y-2)``. The following commands can update
the scripts to the year 2025 using the command line. In future years, increment the
numbers.

.. code-block:: bash

  cd scripts
  sed "s/2024/2025/g" copyrightPrepender.py > tmp
  mv tmp copyrightPrepender.py

  sed "s/2024/2025/g" update_copyright_date.sh > tmp
  sed "s/2023/2024/g" tmp > update_copyright_date.sh
  rm -f tmp
  cd ..


After updating the scripts, commit the changes to a branch since actually running
the update script will cause many files to change.

.. code-block:: bash

  ./scripts/update_copyright_date.sh


Now that the copyrights have been changed, commit the changed files to a branch and
run the next command to prepend a copyright header to files that lack one.

.. code-block:: bash

  python3 scripts/copyrightPrepender.py -r src


Again, commit any changed files to the branch. That is all. Create a pull request
to merge the updated files to the ``develop`` branch.


===================
Web Documentation
===================

Axom web-based documentation is hosted on our 
`Read the Docs project <https://readthedocs.org/projects/axom/>`_. 
Multiple versions are visible there, including the latest content on the 
develop branch (*latest*) and the main branch (*main*). The documentation 
that appears is automatically re-generated each time a change is pushed to 
a branch in the GitHub repository that is enabled to be displayed on the 
Read the Docs project. If you are modifying Axom documentation, you can enable 
the branch you are working on so that you can see what it looks like as you 
push changes to the branch. If your documentation changes are part of a GitHub
pull request, it is a good idea to enable the documentation for that branch
and put a link to it in the pull request summary. This makes it easy for 
reviewers to check over your changes.

.. note :: When you no longer need the documentation of your branch to be
           visible on Read the Docs (e.g., your pull request is merged), 
           please disable that branch on Read the Docs.


===================
Code Health Tools
===================

This section describes how to run code health tools we use.


Code Coverage
---------------

Setting up and running code coverage analysis...


Static Analysis
---------------

Setting up and running static analysis tools....


Memory Checking
----------------

Setting up and running memory checking tools....
