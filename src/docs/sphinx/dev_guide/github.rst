.. ## Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _github-label:

******************************************************
Git/GitHub: Version Control and Branch Development 
******************************************************

This section provides information about getting started with Git and
GitHub and describes some mechanics of topic branch development
on the Axom project. For most project work, we interact with our Git 
repository via our `GitHub project <https://github.com/LLNL/axom>`_.

If you are new to the Git or want to brush up on its features, there are
several good sources of information available on the web:

  * `Atlassian Git Tutorial <https://www.atlassian.com/git/>`_ has a lot of useful stuff.
  * The `Git Docs <https://git-scm.com/docs/>`_ is a complete reference for Git commands and options. It also provides some *cheat sheets* you can download.
  * `Learn Git Branching  <http://learngitbranching.js.org/>`_ is nice for visual, hands-on learners.
  * The e-book `Pro Git, by Scott Chacon <https://git-scm.com/book/en/v2>`_ is an excellent overview guide to using Git effectively.

=========
SSH Keys
=========

If you have not used GitHub before, you should start by creating and adding your SSH keys to GitHub. 
GitHub provides a good tutorial `here <https://help.github.com/en/enterprise/2.18/user/github/authenticating-to-github/adding-a-new-ssh-key-to-your-github-account>`_.
Performing these two simple steps will make it easier for you to interact with 
our Git repository without having to repeatedly enter login credentials. 


.. _repoclone-label:

=========================================
Cloning the Repo
=========================================

All development work on Axom is performed in a *local workspace copy* of
the Git repository. To make a local workspace copy, you clone the repo into 
a directory that you will work in. This is done by typing::

  $ git clone --recursive git@github.com:LLNL/axom.git

.. note:: You don't need to remember the URL for the Axom repo above. It can be
          found by going to the Axom repo on our GitHub project and
          clicking on the 'Clone or download' button on the upper right hand corner
          above the source.

The '--recursive' argument above is needed to pull in all Git *submodules*
that we use in the project. In particular, you will need the BLT build system, 
which is a Git sub-module in Axom, in your local copy of the repo. In case you 
forget to pass the '--recursive' argument to the 'git clone' command, you can
type the following commands after cloning::

  $ cd axom
  $ git submodule init
  $ git submodule update

Either way, the end result is the same and you are good to go.

=========================================
Git Environment Support
=========================================

After cloning, we recommend that you run the development setup script we 
provide in the top-level Axom directory to ensure that your Git environment 
is configured properly; i.e.,::

  $ cd axom
  $ ./scripts/setup-for-development.sh

This script sets up several things we find useful, such as Git editor, aliases,
client-side hooks, useful tips, etc.

You can also define your own aliases for common git commands to simplify your workflow.
For example, the following sets up an alias for unstaging files::

  $ git config alias.unstage 'reset HEAD--'

Then, the alias can be used as a regular Git command as illustrated below::

  $ git unstage <file>
  
Moreover, you may want to tap in to some of your shell's features to enhance your Git experience. 
Chief among the most notable and widely used features are: 

#. Git Completion, which allows tab-completion of Git commands and branch names.

#. Prompt Customization, which allows modifying your prompt to indicate the
   current branch name, whether there are local changes, etc.   

Git ships with contributed plugins for popular shells. Examples illustrating
how to use these plugins in bash and tcsh/csh are given below. 

Setting up your Bash Environment 
--------------------------------

If you are in Bash, you can set your environment as follows:

#. Get the git-prompt.sh and auto-completion scripts from github ::

      $ wget https://raw.githubusercontent.com/git/git/master/contrib/completion/git-prompt.sh
      $ wget https://raw.githubusercontent.com/git/git/master/contrib/completion/git-completion.bash

#. Optionally, you may want to move the files to another location. Nominally, folks put those as 
   hidden files in their home directory ::

      $ mv git-prompt.sh $HOME/.git-prompt.sh
      $ mv git-completion.bash $HOME/.git-completion.bash

#. Add the following to your `.bashrc` ::
  
      source ~/.git-prompt.sh
      source ~/.git-completion.bash
      export GIT_PS1_SHOWDIRTYSTATE=1
      export GIT_PS1_SHOWSTASHSTATE=1
      export GIT_PS1_SHOWUNTRACKEDFILES=1
  
      ## Set your PS1 variable
      reset=$(tput sgr0)
      bold=$(tput bold)
      export PS1='[\w] \[$bold\]$(__git_ps1 " (%s)")\[$reset\]\n\[$bold\]\u@\h\[$reset\] > '
  
Setting up your tcsh/csh Environment 
--------------------------------------

Likewise, if you are using tcsh/csh, you can do the following:

#. Get the auto-completion scripts from github. Note, git-completion.tcsh makes
   calls to git-completion.bash, so you need to have both ::

      $ wget https://raw.githubusercontent.com/git/git/master/contrib/completion/git-completion.tcsh
      $ wget https://raw.githubusercontent.com/git/git/master/contrib/completion/git-completion.bash
  
#. Optionally, you may want to move the files to another location. Nominally, folks put those as
   hidden files in their home directory ::

      $ mv git-completion.tcsh $HOME/.git-completion.tcsh
      $ mv git-completion.bash $HOME/.git-completion.bash
  
#. Add the following to your `.tcshrc` or `.cshrc` ::

      source ~/.git-completion.tcsh
      
      ## Add alias to get the branch
      alias __git_current_branch 'git rev-parse --abbrev-ref HEAD >& /dev/null && echo "{`git rev-parse --abbrev-ref HEAD`}"'
  
      ## Set your prompt variable for example:
      alias precmd 'set prompt="%n@%m[%c2]`__git_current_branch` "'

.. _topicdev-label:

=======================================
Topic Branch Development
=======================================

It is worth re-emphasizing a fundamental principle of the Gitflow
development model that we described in :ref:`gitflow-label`.

.. important:: **We never work directly on the develop or main branches.
               All development occurs on topic branches.**

When we refer to a *topic branch*, it could be a *feature branch*,
a *bugfix branch*, etc. The basic workflow for performing development
on a topic branch is:

  #. Create a topic branch off the develop branch and push the new branch
     to GitHub.
  #. Make changes and commit them to your branch in your local copy of the 
     repository. Remember to push changes to the GitHub repo
     regularly for backup and so you can easily recover earlier versions of
     your work if you need to do so.
  #. If you are working on your topic branch for a while, it is a good idea
     to keep your topic branch current with the develop branch by merging 
     develop into your topic branch regularly. This will simplify the 
     process of merging your work into the develop branch when you are ready.
  #. When your work is complete (including required testing, documentation, 
     etc.), create a pull request so others on the team
     can review your work. See :ref:`pullrequest-label`.

Here are some details about each of these steps.

Step 1 -- Create a topic branch
--------------------------------

    Most development occurs on a topic branch created off the develop branch.
    Occasions where a branch is created from another branch, such as a
    'hotfix' branch created off main, are described in :ref:`gitflow-label`.
    To create a branch in Git, provide the ``-b`` option to the 
    ``git checkout`` command, followed by the name of your topic branch. 
    A topic branch name should include your username (i.e., login id) and a 
    brief description indicating the purpose of the branch. 
    Typically, we label such branches using "feature", "bugfix", etc. to make 
    it clear what type of work is being performed on a branch. For example,::

      $ git checkout -b feature/<userid>/my-cool-new-feature
      $ git push -u

    You can also attach a GitHub issue number to the branch name if the work
    you will do on the branch is related to a issue. Then, GitHub
    will associate the issue with the commit when you merge your branch to the
    develop branch. For example,::

      $ git checkout -b bugfix/<userid>/issue-atk-<issue #>
      $ git push -u

    Alternatively, if your branch addresses multiple issues, you should 
    add the appropriate issue numbers (e.g., #374) to the messages in 
    your commits that address them.

    In each of these examples, the 'git push -u' command pushes the branch to
    the GitHub server and it will appear in the list of branches you and 
    other developers can see there.

Step 2 -- Do development work
--------------------------------

    After you've created a topic branch and pushed it to GitHub, perform 
    your development work on it; i.e., edit files, add files, etc. 
    Common commands you will use are::

      $ git add <file>
      $ git commit
      $ git push

    The 'add' command adds a file (or files) to be staged for a commit 
    operation. The 'commit' command moves your staged changes to your local 
    copy of the repository. The 'push' command pushes these changes to the 
    topic branch in the Git repo. To push your work, you could also do::

      $ git push origin

    This is equivalent to 'git push' if you specified the '-u' option when you
    originally pushed your topic branch when you created it.

 
    .. important:: 
       You may perform several local commits before you push your work 
       to the GitHub repo. Generally, it is a good idea to limit the 
       amount of modifications contained in any one commit. By 
       restricting individual commits to a reasonable size that 
       contain closely related work, it is easier to refer back to 
       specific changes you make when the need arises (as it 
       inevitably will!). For example, if you regularly run your
       code through a formatting tool (we use *clang-format* on the Axom 
       project), it is preferable to commit other content changes first
       and then commit formatting changes in a separate commit. That 
       way, you can distinguish substance from cosmetic changes easily 
       in the Git history.

    Recall the Git environment setup script we recommended that you run after
    cloning the repo in the :ref:`repoclone-label` section above. One of the
    Git pre-commit hooks that the script sets up applies formatting constraints
    on the commit message you provide when you execute the 'commit' command. The
    constraints are recommended Git practices that help make it easier to use
    various tools with the Git version control system. Specifically:

    * Commit message subject line is at most 50 characters
    * Subject line and body of commit message are separated by a blank line
    * Main body of commit message is wrapped to 78 characters

.. _keepcurrent-label:

Step 3 -- Keep current with develop
-------------------------------------

    If you will be working on your topic branch for a while, it is a good idea 
    to merge changes (made by other developers) from the develop branch to 
    your topic branch regularly. This will help avoid getting too far out of 
    sync with the branch into which your work will be merged eventually. 
    Otherwise, you may have many conflicts to resolve when you are ready to 
    merge your topic branch into develop and the merge could be difficult.

    Before you begin the merge, make sure all outstanding changes to your topic
    branch are committed. Then, make sure your local repo is up-to-date with 
    the develop branch by checking it out and pulling in the latest 
    changes; i.e.,::

      $ git checkout develop
      $ git pull

    Next, checkout your topic branch and merge changes in from the
    develop branch, and check for conflicts::

      $ git checkout <your topic branch>
      $ git merge develop

    The 'merge' command will tell you whether there are conflicts and which
    files have them. Hopefully, you will not see any conflicts and you can
    continue working on your topic branch. If there are conflicts, you must
    resolve them before you will be able to merge your topic branch to develop.
    So, you may as well resolve them right away. You can resolve them by
    editing the conflicting files and committing the changes. Merge conflicts
    appear in a file surrounded by lines with special characters on them. For
    example, if you open a conflicted file in an editor, you may see::

      <<<<<<< HEAD
      // lines of code, etc...
      =======
      // more lines of code, etc...
      >>>>>>> develop

    The section above the '=======' line are the file contents in the current 
    branch head (your topic branch). The lines below are the contents of the 
    develop branch that conflict with yours. To resolve the conflict, choose 
    the correct version of contents you want and delete the other lines.

    Alternatively, you can use a tool to help resolve your conflicts. The
    'git mergetool' command helps you run a merge tool. One such tool is called
    "meld", which is very powerful and intuitive. Diff tools like "tkdiff"
    are also helpful for resolving merge conflicts.

    .. important:: **Git will not let you commit a file with merge conflicts.**
                   After you resolve merge conflicts in a file, you must 
                   stage the file for commit (i.e., `git add <filename>),
                   commit it (i.e., `git commit`), and push it to the GitHub
                   repo (i.e., `git push`) before you can merge.

.. _createpr-label:

Step 4 -- Create a pull request
-------------------------------------

    When your work is complete, and you are ready to merge your topic branch 
    to the develop branch, you must initiate a pull request in GitHub. Go
    into the Axom GitHub project, select your branch, and click 
    `Create pull request` in the left column. Make sure you select the correct 
    destination branch. The default destination branch in our project is set 
    up to be the develop branch. So, in most cases, you won't have to do 
    anything special.

    You must also select appropriate team members to review changes. Our 
    GitHub project is set up to require at least one other developer to 
    approve the pull request before a merge.

    .. important:: **You cannot approve your own pull request.**

    When your pull request is approved (see :ref:`review-label` for more
    information), you merge your topic branch to the develop branch by 
    clicking the "merge" button in GitHub. If there are no merge conflicts, 
    the merge will proceed and you are done. If there are conflicts, GitHub 
    will indicate this and will not let you merge until all conflicts are
    resolved.

    The preferred way to resolve conflicts at this point is to go into your 
    topic branch and do the following::

      $ git fetch origin
      $ git merge origin

    The 'fetch' command pulls changes from the remote branch into your local
    branch. Running the 'merge' command will show which files have conflicts.
    Fix the conflicts as described in :ref:`keepcurrent-label`. After all 
    conflicts are resolved, run the 'commit' and 'push' commands as usual::

      $ git commit
      $ git push

    Lastly, complete the merge in GitHub by clicking the merge button.

    .. important:: **To keep things tidy, please delete your topic branch in
                   GitHub after it is merged if you no longer need it for
                   further development. GitHub provides an option to 
                   automatically delete the source branch of a merge after 
                   the merge is complete. Alternatively, you can click on
                   the GitHub branches tab and manually delete the branch.**

================================
Checking Out an Existing Branch
================================

When working on multiple branches, or working on one with someone else on
the team, you will need to checkout a specific branch. Any existing branch
can be checked out from the Git repository. Here are some useful commands::

  $ git fetch
  $ git branch -a
  $ git checkout <branch name>

The 'fetch' command retrieves new work committed by others on branches you may
have checked out, but *without merging* those changes into your local
copies of those branches. You will need to merge branches if you want changes
from one branch to be moved into another. The 'branch' command lists all 
available remote branches. The 'checkout' command checks out the specified 
branch into your local working space.

.. note:: **You do not give the '-b' option when checking out an existing 
          branch. This option is only used when creating a new branch.**

Here is a concrete example::

  $ git branch -a | grep homer
    remotes/origin/feature/homer/pick-up-bart
  $ git checkout feature/homer/pick-up-bart
    Branch feature/homer/pick-up-bart set up to track remote branch feature/homer/pick-up-bart
    Switched to a new branch 'feature/homer/pick-up-bart'

