# Sample .bashrc for SuSE Linux
# Copyright (c) SuSE GmbH Nuernberg

# There are 3 different types of shells in bash: the login shell, normal shell
# and interactive shell. Login shells read ~/.profile and interactive shells
# read ~/.bashrc; in our setup, /etc/profile sources ~/.bashrc - thus all
# settings made here will also take effect in a login shell.
#
# NOTE: It is recommended to make language settings in ~/.profile rather than
# here, since multilingual X sessions would not work properly if LANG is over-
# ridden in every subshell.

# Some applications read the EDITOR variable to determine your favourite text
# editor. So uncomment the line below and enter the editor of your choice :-)
#export EDITOR=/usr/bin/vim
#export EDITOR=/usr/bin/mcedit

# For some news readers it makes sense to specify the NEWSSERVER variable here
#export NEWSSERVER=your.news.server

# If you want to use a Palm device with Linux, uncomment the two lines below.
# For some (older) Palm Pilots, you might need to set a lower baud rate
# e.g. 57600 or 38400; lowest is 9600 (very slow!)
#
#export PILOTPORT=/dev/pilot
#export PILOTRATE=115200


#[[ -f /opt/intel/composer_xe_2013_sp1.2.144/bin/compilervars.sh ]] && . /opt/intel/composer_xe_2013_sp1.2.144/bin/compilervars.sh intel64

[[ -f /opt/intel/composer_xe_2013/bin/compilervars.sh ]] && . /opt/intel/composer_xe_2013/bin/compilervars.sh intel64

[[ -f /opt/intel/impi/4.0.1/bin64/mpivars.sh ]] && . /opt/intel/impi/4.0.1/bin64/mpivars.sh



test -s ~/.alias && . ~/.alias || true


# sasha aliases

alias Fcomp='/data4t/avolkov/Fortran/Comands/comp'
alias Fcompl='/data4t/avolkov/Fortran/Comands/complight'
alias Fcompc='/data4t/avolkov/Fortran/Comands/compc'
alias Fcompo='/data4t/avolkov/Fortran/Comands/compo'

alias Frun='/data4t/avolkov/Fortran/Comands/run'
alias Fmpd='/data4t/avolkov/Fortran/Comands/mpdscr'

# some more ls aliases
alias ll='ls -alFh'
alias la='ls -A'
alias l='ls -CF'
alias lah='ls -lah'
alias Water='cd /data4t/avolkov/Fortran/Shallow_Water/test'
