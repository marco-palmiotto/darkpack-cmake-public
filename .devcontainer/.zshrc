# If not running interactively, don't do anything
[[ $- != *i* ]] && return

# Enable colors and change prompt:
autoload -U colors && colors

# custom prompt
setopt PROMPT_SUBST
PROMPT='%F{magenta}$(date "+%H:%M:%S") %F{red}%n@%M %F{green}%~
%F{blue}%#%f '
TMOUT=1
TRAPALRM() {
    zle reset-prompt
}

# History in cache directory:
HISTSIZE=10000
SAVEHIST=10000
HISTFILE=~/.histfile

# Basic auto/tab complete:
autoload -Uz compinit && compinit
zstyle ':completion:*' menu select
zmodload zsh/complist
compinit
_comp_options+=(globdots)		# Include hidden files.

bindkey -s '^o' 'lfcd\n'

# Edit line in vim with ctrl-e:
autoload edit-command-line; zle -N edit-command-line

# Load aliases and shortcuts if existent.
[ -f "$HOME/.config/shortcutrc" ] && source "$HOME/.config/shortcutrc"
[ -f "$HOME/.config/aliasrc" ] && source "$HOME/.config/aliasrc"

export MAKEFLAGS="-j$(nproc)"

# Load zsh-syntax-highlighting; should be last.
# needs package zsh-syntax-highlighting installed
source /usr/share/zsh/plugins/zsh-syntax-highlighting/zsh-syntax-highlighting.zsh 2>/dev/null



