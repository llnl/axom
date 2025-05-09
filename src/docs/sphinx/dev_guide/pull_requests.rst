.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _pullrequest-label:


********************************
Pull Requests and Code Reviews
********************************

Before any code is merged into the develop or main branches, it
must be tested, documented, reviewed, and accepted. Creating a pull request on
the Axom GitHub project to merge a branch into develop or main 
initiates the test and review processes. All required build configurations
and tests must pass for a pull request to be approved. Also, new tests 
(unit, integration, etc.) must be created to exercise any new functionality 
that is introduced. This will be assessed by reviewers of each pull request. 
See :ref:`createpr-label` for details about creating pull requests.

Code changes in a pull request must be accepted by at least one member
of the Axom development team other than the originator of the pull
request. It is recommended that several team members review pull 
requests, especially when changes affect APIs, dependencies (within Axom
and external), etc. Pull request reviewers can be 
selected on GitHub when the pull request is created. Changes reviewed by 
the team are accepted, rejected, or commented on for improvement; e.g., 
issues to be addressed, suggested changes, etc. Pull requests can be updated
with additional changes and commits as needed. When a pull request is 
approved, it can be merged. If the merged branch is no longer needed for 
development, it should be deleted.

In addition to successful compilation and test passing, changes to the 
develop and main branches should be scrutinized in other ways and using 
other code health tools we use. See :ref:`github-label` for more information 
about using our continuous integration tools.

=======================
Pull Request Summary
=======================

To recap, here is a summary of steps in a pull request:

  #. When code is ready to be considered for acceptance, create a pull request
     on the Axom GitHub project. Identify the appropriate reviewers 
     and add them to the pull request.

  #. Code must build successfully and all relevant tests must pass, including
     new tests required for new functionality.

  #. All issues (build failures, test failures, reviewer requests) must be 
     addressed before a pull request is accepted.

  #. Pull requests must be approved by at least one member of development 
     team other than the pull request originator.

  #. When a pull request is approved it may be merged. If the merged branch is
     no longer needed, it should be deleted. This can be done when merging
     with GitHub. 

.. _review-label:

======================
Code Review Checklist
======================

Beyond build and test correctness, we also want to ensure that code follows
common conventions before acceptance. The following list is a high-level 
summary of the types of concerns we want to identify during pull request 
reviews and resolve before a pull request is merged. Please see the 
`Axom Coding Guide <../coding_guide/index.html>`_ for details
on items in this list.

 #. A new file or directory must be placed in its proper location; e.g.,
    in the same directory with existing files supporting related functionality.
 #. File contents must be organized clearly and structure must be consistent 
    with conventions. 
 #. Namespace and other scoping conventions must be followed. 
 #. Names (files, types, methods, variables, etc.) must be clear, easily
    understood by others, and consistent with usage in other parts of the code.
    Terminology must be constrained; i.e., don't introduce a new term for 
    something that already exists and don't use the same term for different 
    concepts.
 #. Documentation must be clear and follow conventions. Minimal, but adequate, 
    documentation is preferred.
 #. Implementations must be correct, robust, portable, and understandable to
    other developers.
 #. Adequate tests (unit and performance) tests must be added for new 
    functionality.

