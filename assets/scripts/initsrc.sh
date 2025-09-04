
 magic #!, this script must be sourced!

#
#  Begin boilerplate (adapted from setup for shell independence).
#

# Note: All the following special tricks for $_ must continue
#       relaying the value to the next rule.  Be careful!
# Special trick to nail the value of $_ down in a variety of shells.
echo $_ >& /dev/null
# Special trick for tcsh which is one-off on the command history stack.
: $_
# Special trick to capture the value of $_ in zsh and bash
test $?shell$_ != 1$_ >& /dev/null && \
    dollar_underscore="$_" && \
    dollar_underscore=`expr "${dollar_underscore}" : ".\(.*\)"`
# Special trick to capture the value of $_ in tcsh
test $?shell = 1 && set dollar_underscore=`echo $_`

# need to be able to check for mrb
test $?shell = 1 && set ss="csh" || ss="sh"
test "$ss" = "csh" && alias return exit

test "$ss" = "csh" && \
    alias tnotnull "eval '"'test $?'"\!* -eq 1' && eval '"'test -n "$'"\!*"'"'"'"
test "$ss" = "sh" && \
    eval 'tnotnull() { eval "test -n \"\${$1-}\"" ;}'

# Special tricks to figure out if this script has been sourced.
# Works for bash, tcsh, and in some cases for zsh.
set_ is_sourced=false
ifcsh_
    # Note: It is unfortunate that we must hard-code the name
    #       of this script here, but there is no other way
    #       that works, tcsh is brain-dead.
    set base=`basename "${dollar_zed}"`
    test "${base}" != "initsrc.sh" && \
        set is_sourced=true
else
    # Special trick for zsh.
    test "${ZSH_NAME}" && test "${dollar_underscore}" = "${dollar_zed}" && \
        is_sourced=true
    # If there were arguments then there is no safe way to find out
    # whether or not the script was sourced in zsh.  Pretend it was.
    test "${ZSH_NAME}" && test "${#argv}" != "0" && \
        is_sourced=true
    # Special trick for bash.
    test "${BASH}" && test "${BASH_SOURCE}" != "${dollar_zed}" && \
        is_sourced=true
# Warning, this must be here because the tcsh parser is brain-dead.
endif
endifcsh_

#
#  End of boilerplate.  Begin of real work.
#

set_ msg1='ERROR: Could not locate assets dir. Set ASSETS_BASE_DIR.'

# Grid job detection: Set ASSETS_BASE_DIR if in Condor environment with combined tar structure
test "$ss" = "sh" && \
    eval 'if [ -n "${CONDOR_DIR_INPUT-}" ] && [ -d "$CONDOR_DIR_INPUT/strangeness/assets" ]; then
      ASSETS_BASE_DIR="$CONDOR_DIR_INPUT/strangeness/assets"
    fi'
test "$ss" = "csh" && \
    eval 'if ( $?CONDOR_DIR_INPUT && -d "$CONDOR_DIR_INPUT/strangeness/assets" ) then
      setenv ASSETS_BASE_DIR "$CONDOR_DIR_INPUT/strangeness/assets"
    endif'

tnotnull ASSETS_BASE_DIR || ( echo "" ; echo "${msg1}" ; echo "" )
tnotnull ASSETS_BASE_DIR || unset me db dollar_underscore dollar_zed is_sourced base msg1
tnotnull ASSETS_BASE_DIR || return 1

# Original auto-detection logic, adapted for shell independence
if "$ss" = "sh"; then
  if [ -z "${ASSETS_BASE_DIR}" ]; then
    if [ -d "$PWD/assets" ]; then
      ASSETS_BASE_DIR="$(cd "$PWD/assets" && pwd)"
    else
      THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
      CAND="$(cd "$THIS_DIR/.." && pwd)"
      if [ -d "$CAND/calib" ] && [ -d "$CAND/weights" ]; then
        ASSETS_BASE_DIR="$CAND"
      else
        echo "ERROR: Could not locate assets dir. Set ASSETS_BASE_DIR." >&2
        return 1
      fi
    fi
  fi
fi

if "$ss" = "csh"; then
  if ( ! $?ASSETS_BASE_DIR || "$ASSETS_BASE_DIR" == "" ) then
    if ( -d "$PWD/assets" ) then
      set ASSETS_BASE_DIR=`(cd "$PWD/assets" && pwd)`
    else
      set THIS_DIR=`(cd `dirname "$dollar_zed"` && pwd)`
      set CAND=`(cd "$THIS_DIR/.." && pwd)`
      if ( -d "$CAND/calib" && -d "$CAND/weights" ) then
        set ASSETS_BASE_DIR="$CAND"
      else
        echo "ERROR: Could not locate assets dir. Set ASSETS_BASE_DIR." >&2
        return 1
      endif
    endif
  endif
fi

setenv ASSETS_BASE_DIR "${ASSETS_BASE_DIR}"

setenv WEIGHTS_BASE_DIR "${WEIGHTS_BASE_DIR:-$ASSETS_BASE_DIR/weights}"
setenv IA_BADCHANNELS "${IA_BADCHANNELS:-$ASSETS_BASE_DIR/calib/badchannels.txt}"
setenv IA_INFERENCE_WRAPPER "${IA_INFERENCE_WRAPPER:-$ASSETS_BASE_DIR/scripts/run_strangeness_inference.sh}"

setenv PYTHONPATH "${PYTHONPATH:-}:$ASSETS_BASE_DIR:$ASSETS_BASE_DIR/models:$ASSETS_BASE_DIR/scripts"

# Apptainer bindpath handling
if "$ss" = "sh"; then
  if [ -n "${APPTAINER_BINDPATH:-}" ]; then
    export APPTAINER_BINDPATH="$ASSETS_BASE_DIR,${APPTAINER_BINDPATH}"
  else
    export APPTAINER_BINDPATH="$ASSETS_BASE_DIR"
  fi
fi

if "$ss" = "csh"; then
  if ( $?APPTAINER_BINDPATH ) then
    setenv APPTAINER_BINDPATH "$ASSETS_BASE_DIR,${APPTAINER_BINDPATH}"
  else
    setenv APPTAINER_BINDPATH "$ASSETS_BASE_DIR"
  endif
fi

setenv FW_SEARCH_PATH "$ASSETS_BASE_DIR:$ASSETS_BASE_DIR/scripts${FW_SEARCH_PATH:+:$FW_SEARCH_PATH}"

# Report the environment (like setup)
echo
echo ASSETS_BASE_DIR=$ASSETS_BASE_DIR
echo WEIGHTS_BASE_DIR=$WEIGHTS_BASE_DIR
echo IA_INFERENCE_WRAPPER=$IA_INFERENCE_WRAPPER
echo IA_BADCHANNELS=$IA_BADCHANNELS
echo PYTHONPATH=$PYTHONPATH
echo APPTAINER_BINDPATH=$APPTAINER_BINDPATH
echo FW_SEARCH_PATH=$FW_SEARCH_PATH
echo

# Unset temps (like setup)
unset db dollar_underscore dollar_zed is_sourced base msg1
