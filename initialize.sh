#!/usr/bin/env bash

################################################################################

#Setup variables for the script:
UNALIASED_SCRIPT_NAME=$( python -c "import os;print (os.path.realpath(\"${BASH_SOURCE[0]}\"))" )
SCRIPTDIR="$( cd "$( dirname "${UNALIASED_SCRIPT_NAME}" )" && pwd )"
SCRIPTNAME=$( echo $0 | sed 's#.*/##g' )
MINARGS=4
MAXARGS=4
PREREQUISITES="git"

# Determine if this shell is interactive:
ISCALLEDBYUSER=true
[[ "${BASH_SOURCE[0]}" != "${0}" ]] && ISCALLEDBYUSER=false

# Required for any aliased calls:
shopt -s expand_aliases

################################################################################

function simpleUsage()
{
  echo -e "Usage: $SCRIPTNAME PROJECT_NAME \"SHORT PROJECT DESCRIPTION\" \"AUTHOR NAME\" AUTHOR_EMAIL"
  echo -e "Initialize template project variables with project-specific values."
}

#Define a usage function:
function usage() 
{
  simpleUsage
  echo -e ""
  echo -e 'NOTE: This script does not validate the values given, so BE CAREFUL!'
  echo -e ""
  echo -e "Fields that contain spaces must be given in quotes."
  echo -e ""
  # List prerequisites, if any:
  if [[ ${#PREREQUISITES} -ne 0 ]] ; then
    echo -e "Requires the following programs to be installed:"
    for prereq in ${PREREQUISITES} ; do 
      echo -e "  ${prereq}"
    done
    echo
  fi
  echo -e "Return values:"
  echo -e "  0  NORMAL"
  echo -e "  1  TOO MANY ARGUMENTS"
  echo -e "  2  TOO FEW ARGUMENTS"
  echo -e "  3  UNKNOWN ARGUMENT"
  echo -e "  4  MISSING PREREQUISITE"
  echo -e ""
}

#Display a message to std error:
function error() 
{
  echo "$1" 1>&2 
}

# Make a temporary file that will be cleaned after this script finishes. 
function makeTemp()
{
  local f
  f=$( mktemp )
  TMPFILELIST="${TMPFILELIST} $f"
  echo $f
}
TMPFILELIST=''

# Clean all temporary files made in this script:
function cleanTempVars()
{
  rm -f ${TMPFILELIST}
}

# Function run in the exit trap.
function at_exit()
{
  cleanTempVars
}

# Checks the bash built-in PIPESTATUS variable for any failures
# If given strings will output the string corresponding to the failure
# position in PIPESTATUS of any failures in the chain of commands 
# that occurred.
# This function should be called as `checkPipeStatus` as per the 
# alias below it.
function _checkPipeStatus() {

  local hadFailure=false

  for (( i = 0 ; i < ${#RC[@]} ; i++ )) ; do 
    st=${RC[i]}
    echo "st = ${st}"
    if [ $st -ne 0 ] ; then
      # If we were passed a description string for this error in the pipe, then we print it:
      let argIndex=$i+1
      description=${!argIndex}
      [[ ${#description} -ne 0 ]] && error "$description"
      hadFailure=true  
    fi
  done

  if $hadFailure ; then
    return 10
  fi
  return 0
}
alias checkPipeStatus='RC=( "${PIPESTATUS[@]}" );_checkPipeStatus'

function checkPrerequisites {
  local missing=""
  local foundMissing=false
  for c in ${@} ; do
    which $c &> /dev/null
    r=$?
    [[ $r -ne 0 ]] && foundMissing=true && missing="${missing}$c "
  done

  if $foundMissing ; then
    error "Error: the following commands could not be found:"
    error "  $missing"
    error "Please install them and try again."
    if [[ $(uname) == 'Darwin' ]] ; then
      error "NOTE: You likely must use homebrew to install them."
    fi
    return 1
  fi
  return 0
}

# Get an answer from the user in as annoying a method as possible.
# Example:
#  getAnswerFromUser "Are you sure you want to CLEAR the listed records? " 'Yes No' answer 
function getAnswerFromUser()
{
  local prompt="${1}" 
  local acceptableValues="${2}"
  local responseVar="${3}" 

  local haveGoodValue=false
  while ! ${haveGoodValue} ; do
    read -p "${prompt} [$( echo ${acceptableValues} | tr ' ' '/' )] : " ${responseVar}

    for okVal in ${acceptableValues} ; do
      if [[ "$( echo ${!responseVar} | tr a-z A-Z)" == "$( echo ${okVal} | tr a-z A-Z)" ]] ; then
        haveGoodValue=true
      fi
    done

    ! ${haveGoodValue} && error "Please enter one of the following: $( echo ${acceptableValues} | tr ' ' '/' )" && error ""
  done
}


################################################################################

# Do all the interactive stuff if the script is called by a user.
# The bash equivalent to `if __name__ == '__main__':
if ${ISCALLEDBYUSER} ; then
  
  #========================================

  trap at_exit EXIT 
  
  #========================================

  #Check given arguments:
  if [[ $# -gt $MAXARGS ]] ; then
    usage
    exit 1
  elif [[ $# -lt $MINARGS ]] ; then
    usage
    exit 2
  fi

  # Make sure we have all the required commands:
  checkPrerequisites $PREREQUISITES 
  r=$?
  [[ $r -ne 0 ]] && exit 4

  #----------------------------------
  #Read args:

	# Get the variables passed in from the command line.
	# NOTE we do this here because these are all positional arguments
	# and because of the shift operation below.
	proj_name="${1}"
	proj_desc="${2}"
	author_name="${3}"
	author_email="${4}"

  while [ $# -gt 0 ] ; do

    case "$1" in
      -h|--h|--help|-help|help)
        usage;
        exit 0;
        ;;
    esac

    #Get next argument in $1:
    shift
  done

  #----------------------------------
  # Do real work here.

	echo "Populating / renaming files..."

	# Replace trivial fields first.  They should only be in the setup.py:
	tmp_setup=$( makeTemp )
	sed \
		-e "s#PROJECT_NAME#${proj_name}#g" \
		-e "s#_SHORT_PROJECT_DESCRIPTION_#${proj_desc}#g" \
		-e "s#_AUTHOR_EMAIL_#${author_email}#g" \
		-e "s#_AUTHOR_#${author_name}#g" ${SCRIPTDIR}/setup.py > $tmp_setup
	mv $tmp_setup ${SCRIPTDIR}/setup.py

	# Replace fields in the README:
	tmp_readme=$( makeTemp )
	sed \
		-e "s@^# PROJECT\\\\_NAME@# ${proj_name}@" \
		-e "s@^\\\\_SHORT\\\\_PROJECT\\\\_DESCRIPTION\\\\_@${proj_desc}@" ${SCRIPTDIR}/README.md > $tmp_readme
	mv $tmp_readme ${SCRIPTDIR}/README.md

	# Replace project name in tox.ini:
	tmp_tox=$( makeTemp )
	sed \
		-e "s@PROJECT_NAME@${proj_name}@" ${SCRIPTDIR}/tox.ini > $tmp_tox
	mv $tmp_tox ${SCRIPTDIR}/tox.ini

	# Rename the project name in the logger:
	tmp_log=$( makeTemp )
	sed \
		-e "s@PROJECT_NAME@${proj_name}@" ${SCRIPTDIR}/src/PROJECT_NAME/log.py > $tmp_log
	mv $tmp_log ${SCRIPTDIR}/src/PROJECT_NAME/log.py
	
	# Rename the project name in the main executable:
	tmp_main=$( makeTemp )
	sed \
		-e "s@PROJECT_NAME@${proj_name}@" ${SCRIPTDIR}/src/PROJECT_NAME/__main__.py > $tmp_main
	mv $tmp_main ${SCRIPTDIR}/src/PROJECT_NAME/__main__.py

	# Rename the folder containing the code:
	git mv ${SCRIPTDIR}/src/PROJECT_NAME "${SCRIPTDIR}/src/${proj_name}"

	# We're done here:
	echo "Done."
  exit 0
fi

