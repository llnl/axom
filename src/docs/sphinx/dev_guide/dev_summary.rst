.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

****************************************
Development Process Summary
****************************************

This section provides a high-level overview of key Axom software development
topics and includes links to more detailed discussion.


======================================================
Software Development Cycles
======================================================

The Axom team uses a sprint-based development process. We collect
and track issues (bugs, feature requests, tasks, etc.) using ``GitHub``
and define a set of development tasks (i.e., issues) to complete for each 
sprint. While the team meets to discuss issues and plan which ones will be 
worked in each sprint, developers of individual Axom components may plan and 
schedule work in any way that works for them as long as this is coordinated
with other team efforts. Work performed in each sprint work period is tracked 
as a single unified sprint encompassing activities for the entire project.



======================================================
Software Releases and Version Numbers
======================================================

Typically, Axom releases are done when it makes sense to make new features
or other changes available to users. A release may coincide with the completion
of a sprint cycle or it may not.

See :ref:`release-label` for a description of the Axom release process.

The Axom team follows the **semantic versioning** scheme for assigning
release numbers. Semantic versioning conveys specific meaning about 
the code and modifications from version to version by the way version
numbers are constructed.

See :ref:`semver-label` for a description of semantic versioning.


======================================================
Branch Development
======================================================

The Axom project has a ``GitHub`` project space and the team follows 
the **Gitflow** branching model for software development and reviews. Gitflow 
is a common workflow centered around software releases. It makes clear which 
branches correspond to which phases of development and those phases are 
represented explicitly in the structure of the source code repository. As 
in other branching models, developers develop code locally and push their 
work to a central repository.

See :ref:`gitflow-label` for a detailed description of how we use Gitflow.


======================================================
Code Reviews and Acceptance
======================================================

Before any code is merged into one of our main Gitflow branches (i.e., develop 
or main), it must be adequately tested, documented, and reviewed 
for acceptance by other team members. The review process is initiated via 
a *pull request* on the Axom GitHub project.

See :ref:`pullrequest-label` for a description of our review process and 
how we use pull requests.


======================================================
Contributors and Project Access
======================================================

Axom maintains three levels of project access on it GitHub project:

  * **Core team members.** Individuals on the core Axom team are frequent
    Axom contributors and participate regularly in project meetings,
    discussions, and other project activities. They are members of
    the LLNL GitHub organization and the ``axom`` GitHub team. Their
    project privileges include the ability to create branches in the repository,
    push code changes to the Axom repo, make PRs, and merge them when they are
    approved and all checks have passed.
  * **Regular contributors.** Individuals, who are not on the core Axom team,
    but are members of the LLNL GitHub organization and are involved in some
    aspects of Axom development are considered regular contributors. They are
    members of the ``axom-contrib`` GitHub team. Their project privileges
    include the ability to create branches in the repository, push code changes
    to the Axom repo, and make PRs. However, they may not merge PRs and must
    coordinate with the core team to have their work included in the develop
    branch. This is mainly due to the way GitHub structures its project
    access levels.
  * **Everyone else.** Anyone with a GitHub account is welcome to contribute
    to Axom. Individuals outside of the two groups described above, and 
    specifically not a member of LLNL GitHub organization, can make PRs
    in the Axom project, but must do so from a branch on a *fork* of
    the Axom repo. Thus, the process of reviewing and merging contributions
    involves additional steps which we describe here.

--------------------------
Forking the repository
--------------------------

The requirement for individuals outside of the LLNL GitHub organization
to contribute on a fork of the repo is due to policies enforced
by the LLNL organization on GitHub (in which the Axom project resides) and the
Livermore Computing (LC) organization (in which we run our GitLab CI testing).
Fortunately, you may still contribute to Axom by `forking the Axom repo
<https://github.com/LLNL/axom/fork>`_. Forking creates a copy of the Axom
repository that you own. You can make changes on your local copy and push them
to your fork on GitHub. When you are ready to have your Axom contribution
reviewed and added to the Axom project, you may create a pull request in the 
Axom project.

--------------------------------------------------
Accepting a pull request from a forked repository
--------------------------------------------------

Due to LLNL security policies, some Axom pull requests cannot be run through 
all Axom CI checks. The Livermore Computing (LC) Center GitLab systems 
restrict which GitHub PRs may run automatically through its CI test pipelines.
For example, a PR made from branch on a forked repository will not trigger 
GitLab CI checks. GitLab CI on LC platforms will be run only on PRs that are 
made from branches in the GitHub Axom repository.

.. note:: **The following process for accepting PR contributions from a fork
          of the Axom repo must be executed by a member of the Axom team:**

          To facilitate testing contributions in PRs from forked repositories,
          we maintain a script to pull a PR branch from a forked repo into the
          Axom repo. First, identify the number of the PR, which appears at
          the top of a PR. Then, run a script from the top-level Axom
          directory::

            $ ./scripts/make_local_branch_from_fork_pr.sh -b <PR #>

          If successful, this will create a branch in your local copy of the
          Axom repo labeled ``pr-from-fork/<PR #>`` and you will be on that
          local branch in your checkout space. To verify this, you can run
          the following command after you run the script::

            $ git branch

          You will see the new branch in the listing of branches and the branch
          you are on will be starred.

          You can push the new branch to the Axom repo on GitHub::

            $ git push git@github.com:LLNL/axom.git <branch-name>

          and make a PR for the new branch. It is good practice to reference
          the original PR in the description of the new PR to track the
          original PR discussion and reviews.

          All CI checks will be triggered to run on the new PR made in the
          Axom repo. When everything passes and the PR is approved, it may
          be merged. When it is merged, the original PR from the forked repo
          will be closed and marked as merged unless it is referenced
          elsewhere, such as in a GitHub issue. If this is the case, then the
          original PR (from the forked repo) must be closed manually.


======================================================
Testing and Code Health
======================================================

Comprehensive software testing processes and use of code health tools (e.g., 
static analysis, memory checkers, code coverage) are essential ingredients 
in the Axom development process.

See :ref:`testing-label` for a description of our software testing process,
including *continuous integration*.


======================================================
Software Development Tools
======================================================

In addition to the tools listed above, we use a variety of other tools to help
manage and automate our software development processes. The *tool philosophy*
adopted by the Axom project focuses on three central tenets:

  * Employ robust, commonly-used tools and don't re-invent something that already exists.
  * Apply tools in ways that non-experts find easy to use.
  * Strive for automation and reproducibility.

The main interaction hub for Axom developers is the **Atlassian
tool suite** on the Livermore Computing Collaboration Zone (CZ). These tools
can be accessed through the `MyLC Portal <https://lc.llnl.gov>`_.
Developer-level access to Axom project spaces in these tools requires 
membership in the LC group 'axom'. If you are not in this group, and need 
to be, please send an email request to 'axom-dev@llnl.gov'.

The main tools we use are listed below. Please navigate the links
provided for details about how we use them and helpful information about 
getting started with them.

* **Confluence.**  We use the `Axom Confluence space <https://lc.llnl.gov/confluence/display/ASCT>`_ for team discussion (e.g., hashing out design ideas), maintaining meeting notes, etc.

* **GitHub.** We use the `Axom GitHub project <https://github.com/LLNL/axom>`_ to manage our issues and Git repository which contains the Axom source code, build configurations, scripts, test suites, documentation, etc.

  * See :ref:`github-label` for more information about how we use Git and GitHub.

* **GitLab.** We use GitLab for continuous integration to ensure code quality on our LC systems.:  `Axom GitLab project <https://lc.llnl.gov/gitlab/axom/axom>`_

  * See :ref:`gitlab-label` for more information about how we use GitLab.

* **GitHub Actions.** We use GitHub Actions for continuous integration to ensure every code change passes a
  level of quality before being merged.:  `GitHub Actions <https://github.com/LLNL/axom/actions>`_

  * See :ref:`github_actions-label` for more information about how we use GitHub Actions.
