#!/usr/bin/env bash

# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#------------------------------------------------------------------------------
# This script sets up aliases to simplify development.
#------------------------------------------------------------------------------

read -ep 'Would you like to configure git aliases? [Y/n] ' configure &&
if test "$configure" = "n" -o "$configure" = "N"; then
    exit 0
fi

## shortcut for git config
gconf="git config"

## Usage: git pullall
##
## Pulls latest including the latest submodules
${gconf} alias.pullall "!bash -c \"git pull && git submodule update --init --recursive\""

## Usage: git unstage <file>
## 
## Unstages the given file.
${gconf} alias.unstage "reset HEAD --"

## Usage: git find <string>
##
## Finds all the files matching the given string.
${gconf} alias.find "!bash -c \"git ls-tree -r --name-only HEAD | grep --color \$1\""

## Usage: git [incoming] [outgoing]
## 
## Lists incoming/outgoing commits from the current branch to develop.
${gconf} alias.incoming "!bash -c \"git fetch origin && git log --pretty=format:'%C(bold)%h%Creset %s %C(ul)(%an)%Creset' HEAD..origin/develop\""
${gconf} alias.outgoing "!bash -c \"git fetch origin && git log --pretty=format:'%C(bold)%h%Creset %s %C(ul)(%an)%Creset' origin/develop..HEAD\""

## Usage: git branch-history
##
## Displays the history graph of the current branch
${gconf} alias.history-graph "log --graph --pretty=format:'%C(bold yellow)%h%Creset %ad %C(bold normal)%s%Creset %C(yellow)%d%Creset %C(bold blue)<%an>%Creset' --date=short"

## Usage: git all-branch-history
##
## Displays the history graph of all the branches
${gconf} alias.history-graph-all "log --graph --full-history --all --pretty=format:'%C(bold yellow)%h%Creset %ad %C(bold normal)%s%Creset %C(yellow)%d%Creset %C(bold blue)<%an>%Creset' --date=short"

## Usage: git show-config
##
## Displays the git configuration settings
${gconf} alias.show-config "!bash -c \"git config --list\""

## Usage: git show-aliases
##
## Displays the current set of git aliases.
${gconf} alias.show-aliases "!bash -c \"git config --list | grep alias | sort\""

## Usage:: git prepush
##
## Provides a snapshot of the changes to be pushed
${gconf} alias.prepush "log --graph --stat origin/develop.."

## Usage: git undo
##
## An alias to undo the last `git commit`, but leave
## the unstaged changes in the file-system.
${gconf} alias.undo "!bash -c \"git reset HEAD~ \""

## Usage: git squash-merge <branch-name>
##
## Merges the branch with the given name into the current branch,
## typically main, using a squash-merge. All commits in the
## branch <branch-name> are squashed to a single commit.
${gconf} alias.squash-merge "!bash -c \"git merge --squash \$1\""

## Usage: git ffwd-merge <branch-name>
## 
## Does a fast-forward merge of the given branch into
## the current branch, typically main. All commits in
## branch <branch-name> are appended/played on top of
## main. 
##
## A precondition for a ffwd-merge is that a linear 
## commit history is maintained as depicted in the figure
## below:
##
##			  D---E---F <branch-name>
##			 /
##  A---B---C <main>
##
## In order to maintain the topic branch in this state, 
## `git rebase main` must be issued before the
## `git ffwd-merge <branch-name>`.
${gconf} alias.ffwd-merge "!bash -c \"git merge --ff-only \$1\""
