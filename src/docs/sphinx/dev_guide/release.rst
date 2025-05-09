.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _release-label:

*******************************************
Release Process
*******************************************

The Axom team decides as a group when the code is ready for a release.
Typically, a release is done when we want to make changes available to users;
e.g., when some new functionality is sufficiently complete or we want users to
try something out and give us feedback early in the development process. A
release may also be done when some other development goal is reached. This
section describes how an Axom releases is done. The process is fairly
informal. However, we want to ensure that the software is in a reasonably
robust and stable state when a release is done. We follow this process to
avoid simple oversights and issues that we do not want to pass on to users.

In the :ref:`gitflow-label` section, we noted that the **main branch
records the official release history of the project**. Specifically,
whenever, the main branch is changed, it is tagged with a new
version number. We use a git 'lightweight tag' for this purpose. Such
a tag is essentially a pointer to a specific commit on the main branch.

We finalize preparations for a release on a release branch so that other
work may continue on the develop branch without interruption.

.. note:: No significant code development is performed on a release branch.
          In addition to preparing release notes and other documentation, the
          only code changes that should be done are bug fixes identified
          during release preparations

Here are the steps to follow for an Axom release.

1: Start Release Candidate Branch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a release candidate branch off of the develop branch to initiate a
release. The name of a release branch should contain the associated release
version name. Typically, we use a name like v0.5.0-rc
(i.e., version 0.5.0 release candidate). See :ref:`semver-label` for a
description of how version numbers are chosen.

2: Create a Pull Request
^^^^^^^^^^^^^^^^^^^^^^^^

The release candidate branch, when complete, reviewed, and approved, will be
merged into main so that a release tag can be generated for that merge commit.
Create a pull request to merge the release candidate branch into main so that
release changes can made and reviewed. 

.. note:: Typically, when a release is being prepared it will have been months
          since the previous release and the main branch has changed. Thus,
          the number of files changed by the release candidate merge into main
          will be large. Fortunately, most of those changes have been reviewed
          and merged into the develop branch and do not require
          additional review. Therefore, it is helpful to create a companion
          pull request to merge the release candidate branch into develop. 
          This pull request will not be merged, but will be much easier for
          the team to review. To facilitate the process, cross reference the
          pull request to be merged into main and the one to develop and note 
          in the summary of the former to review the latter but approve the
          former. (Whew! Hopefully, that is clear!)

Typical changes that are made in a release candidate branch include:

#. Update the version information (major, minor, and patch version numbers)
   at the top of the ``axom/src/cmake/AxomVersion.cmake`` file and in
   the ``axom/RELEASE`` file.

#. Update the notes in the section for the new release in the file 
   ``axom/RELEASE-NOTES.md``. Add the release version number and release date
   in the section heading and add a link to the new version on GitHub at the
   bottom of the file.

#. Update the mail map in ``axom/.mailmap``, if needed, by adding names and 
   emails of new contributors since the last release.

#. Update the citations in ``axom/CITATION.cff`` by adding the names
   of new LLNL contributors since the last release.

#. Test the code by running it through all continuous integration tests
   and builds. This will be done automatically when the release pull request is
   made. All build configurations must compile properly and all tests must pass
   before the pull request can be merged.

#. Fix any issues discovered during final release testing in the release 
   candidate branch if code changes are reasonably small, and re-run
   appropriate tests to ensure issues are resolved. If a major bug is
   discovered, and it requires significant code modifications to fix,
   do not fix it on the release branch.  `Create a new GitHub issue for it
   <https://github.com/LLNL/axom/issues/new>`_ and note it in the ``known bugs``
   section of the release notes. Alternatively, if time permits, fix the 
   bug in a different branch and create a pull request as you would do during
   regular development. After the bug is resolved and that pull request is
   merged into develop, merge develop into the release candidate branch where
   checks will run on that.

#. Make sure all documentation (source code, user guides, etc.) is
   updated and reviewed. This should not be a substantial undertaking as
   most of this should have been done during the regular development cycle.

#. Proofread the release notes for completeness and clarity and address
   any shortcomings. Again, this should not take much time as release notes
   should be updated during the regular development cycle. See
   :ref:`release-notes-label` for information about release notes.

.. important:: It is good practice to have everyone on the team review the
               release notes to ensure that they are complete, correct, and
               sufficiently descriptive so that users understand the content
               of the release. **Please make sure the section for the new
               release follows the same organization as in previous release
               sections.**

3: Merge Release Candidate
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Merge the release candidate branch into the main branch once it is ready and
approved. Do not "squash merge" as it will make the histories of main and
develop branches disagree, and we want to preserve the history. After
merging, the release candidate branch can be deleted.


4: Draft a GitHub Release
^^^^^^^^^^^^^^^^^^^^^^^^^

`Draft a new Release on GitHub <https://github.com/LLNL/axom/releases/new>`_

#. Enter the desired tag version, e.g., v0.5.0

#. Select **main** as the target branch to tag a release.

#. Enter a Release title. We typically use titles of the following form *Axom-v0.3.1*

#. Copy and paste the information for the release from the
   ``axom/RELEASE-NOTES.md`` into the release description (omit any sections if empty).

#. Publish the release. This will create a tag at the tip of the main
   branch and add corresponding entry in the
   `Releases section <https://github.com/LLNL/axom/releases>`_

.. note::

   GitHub will add a corresponding tarball and zip archives consisting of the
   source files for each release. However, these files do not include any
   submodules. Consequently, a tarball that includes all of the submodules is
   generated manually in a separate step.

5: Make a Release Tarball
^^^^^^^^^^^^^^^^^^^^^^^^^^

* Checkout the main branch locally and run ``axom/scripts/make_release_tarball.sh --with-data``
  Run this script from the top level ``axom`` subdirectory. This will
  generate a two tarballs of the form ``Axom-v0.3.1.tar.gz`` and ``AxomData-v0.3.1.tar.gz``
  consisting of the axom source and data respectively.

* Upload the tarballs for the corresponding release, by going to the
  `GitHub Releases section <https://github.com/LLNL/axom/releases>`_ and ``Edit``
  the release created earlier.

* Attach the tarball to the release.

* Add a note at the top of the release description that indicates which
  tarball consists of all the submodules, e.g., *\"Please download the Axom-v0.3.1.tar.gz tarball below, which includes all of the Axom submodules as well\"*

* Update the release.

6: Merge Main to Develop
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a pull request to merge main into develop so that changes in the 
release candidate branch are integrated into subsequent Axom development.
When approved, merge it.


.. _release-notes-label:

*******************************************
Release Notes
*******************************************

Axom release notes are maintained in a single file ``axom/RELEASE-NOTES.md``.
The release notes for the latest version are at the top of the file.
Notes for previous releases appear after that in descending version number
order.

For each version, the release notes must contain the following information:

 * Axom version number and date of release

 * One or two sentence overview of release, including any major changes.

 * Release note items should be broken out into the following sections:

    * Added: Descriptions of new features
    * Removed: Notable removed functionality
    * Deprecated: Deprecated features that will be removed in a future release
    * Changed: Enhancements or other changes to existing functionality
    * Fixed: Major bug fixes
    * Known bugs: Existing issues that are important for users to know about

.. note:: Release notes for each Axom version should explain what changed in
          that version of the software -- and nothing else!!

Release notes are an important way to communicate software changes to users
(functionality enhancements, new features, bug fixes, etc.). Arguably, they
are the simplest and easiest way to do so. Each change listed in the release
notes should contain a clear, concise statement of the change. Items should
be ordered based on the impact to users (higher impact - first, lower impact
last).

.. note:: When writing release notes, think about what users need to know and
          what is of value to them.

Release notes should summarize new developments and provide enough detail
for users to get a clear sense of what's new. They should be brief -- don't
make them overly verbose or detailed. Provide enough description for users
to understand a change, but no more than necessary. In other words, release
notes summarize major closed issues in a human-readable narrative. Direct
users to other documentation (user guides, software documentation, example
codes) for details and additional information.

.. note:: Release notes should be updated as work is completed and reviewed
          along with other documentation in a pull request. This is much
          easier than attempting to compile release notes before a release
          by looking at commit logs, etc. Preparing release notes as part
          of the release process should take no more than one hour.

Lastly, release notes provide an easy-to-find retrospective record of
progress for users and other stakeholders. They are useful for developers
and for project reporting and reviews.


