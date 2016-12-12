#!/usr/local/bin/python
#
# This modules contain helper functions to submit SGE blast and blat jobs.
#
# Author: Michelle Dimon

# ORIGINAL IMSA_v2 file
# Protected by license (see IMSA_License.txt)

import subprocess

import config

def runProcess(exe):
    """Runs a command using subprocess.Popen, returns stdout and stderr strings."""
    p = subprocess.Popen(exe, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out,err = p.communicate()

    return out, err

def getQstatMsg():
    """Runs qstat and returns the qstat message string."""
    qstatMsg, qstatErr = runProcess("qstat")

    if len(qstatMsg) < 2:
        print "Error retrieving qstat message, or no running jobs."
        print "msg was:", qstatMsg
        print "err was:", qstatErr
        qstatMsg, qstatErr = runProcess("qstat")

    return qstatMsg

def getRunningJobNames():
    """Uses qstat to determine the names of the jobs that are running."""
    qstatMsg = getQstatMsg()

    names = []
    for line in qstatMsg.split("\n"):
        #print "line:", line
        if line.startswith("job-ID") or line.startswith("--") or len(line) < 2:
            continue
        pieces = line.split()
        names.append(pieces[2])

    return names

def getRunningJobIds():
    """Uses qstat to determine the ids of the jobs that are running."""
    qstatMsg = getQstatMsg()

    ids = []
    for line in qstatMsg.split("\n"):
        if line.startswith("job-ID") or line.startswith("--") or len(line) < 2:
            continue
        pieces = line.split()
        ids.append(pieces[0])

    return ids

def isJobRunning(jobname):
    """Uses qstat to determine whether the given jobname is in the set of jobs
    that are currently running."""
    qstatMsg = getQstatMsg()

    return jobname in qstatMsg


def createBlastScript(scriptFilename, jobname, db, query, output, params=""):
    """Given the blast details (db, query, other params), this function creates as script
    that will run the blast job."""

    if len(jobname) > 10:
        raise Exception("The job name can only be 10 characters or less.  Otherwise it can't be monitored")
    out = open(scriptFilename, "w")
    out.write("""
#! /bin/bash
### Change to the current working directory:
#$ -cwd
### Job name:
#$ -N %s

%s -outfmt 6 %s -db %s -query %s -out %s

""" % (c.PATH_TO_BLASTN, jobname, params, db, query, output))

def startBlastJob(jobname, db, query, output, params):
    """Given the blast details (db, query, other params), this function creates a script
    that will run the given blast job and then uses qsub to submit the job to SGE."""
    scriptname = "run"+jobname+".sh"
    createBlastScript(scriptname, jobname, db, query, output, params)

    cmd = "qsub %s" % (scriptname)
    try:
        p = subprocess.Popen([cmd], shell=True)
        p.wait()
    except:
        raise

def createBlatScript(scriptFilename, jobname, db, query, output, params=""):
    """Given the blat details (db, query, other params), this function creates a script
    that will run the given blat job."""
    if len(jobname) > 10:
        raise Exception("The job name can only be 10 characters or less.  Otherwise it can't be monitored")
    out = open(scriptFilename, "w")
    out.write("""
#! /bin/bash
### Change to the current working directory:
#$ -cwd
### Job name:
#$ -N %s

blat %s %s %s %s

""" % (jobname, db, query, params, output))

def startBlatJob(jobname, db, query, output, params):
    """Given the blat details (db, query, other params), this function creates a script
    that will run the given blat job and then uses qsub to submit the job to SGE."""
    scriptname = "run"+jobname+".sh"
    createBlatScript(scriptname, jobname, db, query, output, params)

    cmd = "qsub %s" % (scriptname)
    try:
        p = subprocess.Popen([cmd], shell=True)
        p.wait()
    except:
        raise



                                    
