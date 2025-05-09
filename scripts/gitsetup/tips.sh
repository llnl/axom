#!/usr/bin/env bash

# Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

##==============================================================
# This script suggests optional settings that developers
# may want to set, such as a git prompt, auto-completion,
# default editor, etc.
##==============================================================

read -ep 'Would you like to see a list of Git tips/suggestions? [Y/n] ' seetips &&
if test "$seetips" = "n" -o "$seetips" = "N"; then
    exit
fi

echo "Here is a list of tips/suggestions that may be helpful when working with Git..."
## Suggest how to set the default editor
echo '
******************************
*       EDITOR SETUP         *
******************************
The default editor used to compose git commit messages may be set with:

git config --global core.editor <your-editor>

For example, for setting up VIM, you would do:

    git config --global core.editor vim

Likewise, for setting up emacs, you would do:

    git config --global core.editor emacs


'

## Suggest setting up a merge-tool
current_merge_tool="$(git config --get merge.tool)"
if test -n "$current_merge_tool"; then
  echo "You are currently using: $current_merge_tool"
fi


echo '
******************************
*       MERGETOOL SETUP      *
******************************
One can configure git to to load a particular merge tool with

    git config --global merge.tool <toolname>


To see a list of available mergetools run the following:

    git mergetool --tool-help

Additional information is available "git help mergetool".


'

## Suggest setting up prompt and command tab-completion
echo '
******************************
* PROMPT/GIT COMPLETION      *
******************************
You may setup your prompt to do tab auto-completion of git commands
by sourcing the appropriate git-completion script for your shell:

    https://github.com/git/git/tree/master/contrib/completion

In addition, you may want to setup your prompt to dynamically show you
on what branch you are on, if you have any stashed/staged modifications or
any untracked files.

To do this you can just get and follow the instructions in the beginning
of the git-prompt.sh script that is available here:

    https://github.com/git/git/blob/master/contrib/completion/git-prompt.sh

This only supports Bash/zsh shell environments.

If you are using c-shell/tcsh you may want to follow the instructions in the
link below on how to setup your git-prompt in tcsh:

    http://thrysoee.dk/gittcsh/


'
