#!/usr/bin/env bash

# ngsadmix.sh - a wrapper to help ensure convergence of ngsadmix replicates.
# Execute with whatever arguments you would give to ngsadmix and it will run
# replicates until they converge or a maximum number is reached

# Depends on: NGSadmix (case sensitive) in $PATH

# Author: Zachary J. Nolen
# Last updated: 2023-04-18


# If you would like to change the settings of this wrapper, export the
# variables of interest before running and it will prioritize the ones you
# define. If you don't, they will set their defaults. Eligible options are:
#
# $TMPDIR - The temporary directory where temporary results are processed
#
# $reps - The maximum number of replicates to perform if convergence is not
# reached (default 100)
#
# $minreps - The minimum number of replicates to perform before assessing if 
# convergence is reached (default 10)

# $conv - The number of top replicates to assess convergence of (default 3)
#
# $thresh - The maximum range of likelihoods for the top $conv replicates to 
# count as converged (default 2)
#
# Defaults are adapted from methods performed in Pečnerová et al. 2021, Current
# Biology. 

#####################
### Housekeeping ####
#####################

# Define function to make printing to stderr easy
echoerr() { echo "$@" 1>&2; }

# Make it so stuff fails more often
set -eo pipefail

#####################
## Parse Arguments ##
#####################

# Parse arguments and divide up between what will be passed to ngsadmix 
# and what will be used for wrapper processing. Method adapted from this
# useful answer on stack overflow: https://stackoverflow.com/a/14203146

passedargs=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -o|-outfiles)
      outpre="$2"
      shift # past argument
      shift # past value
      ;;
	-resume)
	  resume=1
	  shift
	  ;;
    -*|--*)
      passedargs+=("$1 $2")
	  shift
	  shift
      ;;
    *)
      passedargs+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

passedargs=${passedargs[@]}

# log what happened
echoerr "The following options will be passed to ngsadmix, please make sure"
echoerr "they are correct:"
echoerr "$passedargs -o $outpre"
echoerr

#####################
# Setup directories #
#####################

outfile=$(basename $outpre)
outdir=$(dirname $outpre)

# temp directory
if [ -z $TMPDIR ]; then
    TMPDIR=$HOME/scratch
fi
temp=$TMPDIR/ngsadmix_opt_$outfile
mkdir -p $temp

# out directory
mkdir -p $outdir

# Make sure output and temp directories aren't the same. There is a very low 
# chance they would be, but why risk it...
if [ $temp == $outdir ]; then
	echoerr "ERROR: Temporary directory is the same as output directory."
	echoerr 'Please change variable $TMPDIR or -o/-outfiles path so that'
	echoerr "they do not match."
	exit -1
fi

# log what happened
echoerr "Temporary files will be written to:"
echoerr $temp
echoerr
echoerr "Highest likelihood replicate will be written with this prefix:"
echoerr $outpre
echoerr

#####################
##### Set params ####
#####################

# set maximum number of replicates to perform before stopping if not 
# converged
if [ -z $reps ]; then
	reps=100
else
	echoerr 'WARNING: $reps has been modified from the default.'
	echoerr 'Ignore if this change was intentional.'
	echoerr
fi

if [ -z $minreps ]; then
	minreps=20
else
	echoerr 'WARNING: $minreps has been modified from the default.'
	echoerr 'Ignore if this change was intentional.'
	echoerr
fi

# set maximum range between top reps in likelihood units to qualify reps
# as reaching convergance
if [ -z $thresh ]; then
	thresh=2
else
	echoerr 'WARNING: $thresh has been modified from the default.'
	echoerr 'Ignore if this change was intentional.'
	echoerr
fi

# set the minimum number of top replicates that must be within the value 
# $thresh to qualify as converged
if [ -z $conv ]; then
	conv=3
else
	echoerr 'WARNING: $conv has been modified from the default.'
	echoerr 'Ignore if this change was intentional.'
	echoerr
fi

# log what got set
if (( $(echo "$thresh > 0" | bc -l) )); then
echoerr "Beginning ngsadmix replicates, will perform replicate runs until "
echoerr "until $reps replicates have finished. If you'd like to change "
echoerr 'these convergence criteria, edit the $reps, $thresh, or $conv '
echoerr 'variables in the wrapper or export them before starting.'
echoerr
else
echoerr "Beginning ngsadmix replicates, will perform replicate runs until "
echoerr "the top $conv runs are within $thresh log-likelihood units of each "
echoerr "other or until $reps replicates have finished. If you'd like to "
echoerr 'change these convergence criteria, edit the $reps, $thresh, or '
echoerr '$conv variables in the wrapper or export them before starting.'
echoerr
fi

#####################
## Run the wrapper ##
#####################

# set log file where replicate likelihoods will be stored and assessed from
templog=$temp/${outfile}_optimization_wrapper.log
permlog=$outdir/${outfile}_optimization_wrapper.log
mkdir -p $temp/bestrep

# initialize the log and make the temporary directory for kept results
if [ "$resume" == "1" ]; then
	echoerr "Resuming optimization from previous run."
	if [ ! -f $permlog ]; then
		echoerr "ERROR: Log from previous optimization not found, please "
		echoerr "ensure that an optimization log from a previous run of this "
		echoerr "wrapper, along with its best replicate output are in: "
		echoerr "$permlog"
		exit -1
	fi
	cp $permlog $templog
	cp $outdir/${outfile}.filter $temp/bestrep/
	cp $outdir/${outfile}.fopt.gz $temp/bestrep/
	cp $outdir/${outfile}.log $temp/bestrep/
	cp $outdir/${outfile}.qopt $temp/bestrep/
	initrep=$(tail -1 $templog | awk '{print $1+1}')
	echoerr "Will pick up on replicate $initrep, if run has not already "
	echoerr "met convergence criteria."
	echoerr "WARNING: I will not check to ensure your options are the same "
	echoerr "on this run as the previous, only that the output filenames "
	echoerr "match."
else
	> $permlog
	> $templog
	initrep=1
fi

# run replicates of ngsadmix until we hit either convergence or $reps
for i in $(seq $initrep $reps); do
	if (( $(echo "$thresh > 0" | bc -l) )); then
		# check if loop should break when convergence is reached
		if (( "$i" > "$minreps" )); then
			# sort likelihoods, see if top three are within $thresh of each other
			# and break if so
			diff=$(sort -k3gr $templog | head -n 3 | awk '{print $3}' | \
					tr '\n' ' ' | awk '{print $1-$3}')
			if (( $(echo $diff $thresh | awk '{print $1 <= $2}') )); then
				break
			fi
		fi
	fi

	# make a folder for the replicate results
	reptemp=$temp/iter$i
	mkdir -p $reptemp

	# run replicate
	echoerr "Beginning replicate $i..."
	seed=$RANDOM
	NGSadmix $passedargs -o $reptemp/$outfile -seed $seed

	# add likelihood to log
	like=$(grep 'best like=' $reptemp/$outfile.log | sed 's/=/ /g' | \
			awk '{print $3}')
	echo "$i	$seed	$like" >> $templog
	cp $templog $permlog

	# if better than previous reps, store results
	bestrep=$(sort -k3gr $templog | head -n 1 | awk '{print $1}')
	if [ "$bestrep" == "$i" ]; then
		echoerr "Updating best replicate..."
		echoerr "WARNING: If the script is cancelled before I say 'Done!' "
		echoerr "there may be a mismatch in the optimization log and saved "
		echoerr "best replicate."
		mv $reptemp/* $temp/bestrep/
		mv $temp/bestrep/* $outdir/
		echoerr "Done!"
	fi

	# delete replicate folder
	rm -r $reptemp

done

diff=$(sort -k3gr $permlog | head -n 3 | awk '{print $3}' | tr '\n' ' ' | \
	awk '{print $1-$3}')

if (( $(echo $diff $thresh | awk '{print $1 <= $2}') )); then
	echoerr
	echoerr "Optimization complete! Top $conv replicates have converged with"
	echoerr "$diff log-likelihood unit difference."
	echoerr
else
	echoerr
	echoerr "WARNING: Optimization reached max reps ($reps) without converging"
	echoerr "with a difference of $diff log-likelihood units between the top"
	echoerr "$conv results. Consider increasing the number of replicates or"
	echoerr "performing the analysis only for lower values of K."
	echoerr
fi

echoerr "Final checks to make sure you have the right files..."
bestlike=$(sort -k3gr $templog | head -n 1 | awk '{print $3}')
yourlike=$(grep 'best like=' $outpre.log | sed 's/=/ /g' | \
			awk '{print $3}')

if (( $(awk 'BEGIN{ print "'$bestlike'"=="'$yourlike'" }') )); then
	echoerr "Best likelihood results successfully written to output directory!"
else
	echoerr "ERROR: Something went wrong, try comparing the values of the"
	echoerr "optimization_wrapper.log file and the .log file of the ngsadmix"
	echoerr "run."
	echoerr "best like = $bestlike; output like = $yourlike"
	exit -1
fi