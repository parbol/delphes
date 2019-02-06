##############################################################################
##############################################################################
##                                                                          ##
##                              submit.py                                   ##
##                                                                          ##
##            Tool for submitting Delphes jobs to lxbatch.                  ##                        
##    Check: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MTDPhysicsTDR     ##
##############################################################################
##############################################################################
import os, sys, optparse
import ROOT as r

templateLSF = """
#!/bin/bash
pushd CMSSWRELEASE/src
eval `scramv1 runtime -sh`
pushd
EXENAME CARD OUTPUT FIRSTEVENT LASTEVENT INPUT 
"""


templateIFCA = """

###################
# FLAG definition #
###################

# Request the Bourne Shell
#$ -S /bin/bash

# Change to the current working directory before starting the job
#$ -cwd

# Change the job name to "hello_world"

#$ -N NAMEOFJOB 

# Resource request. We request 1MB of memory, and 60 seconds of wall
# clock time, that more than is enough for the test.
# -l mem_free=2G
# -l h_rt=4:0:00

# We are using the "l.tests" project for the examples
#$ -P l.gaes


#################
# Actual script #
#################

source /cvmfs/cms.cern.ch/cmsset_default.sh
/cvmfs/cms.cern.ch/common/scramv1
pushd CMSSWRELEASE/src
eval `scramv1 runtime -sh`
pushd
EXENAME CARD OUTPUT FIRSTEVENT LASTEVENT INPUT 
"""



##############################################################################################################################################
def prepareJob(exeName, fileName, initevent, endevent, outputfilename, logfilename, errfilename, jobname, cardLocation, cmsswrelease, launcher, queue, site):

    thisText = templateLSF
    if site == 'ifca':
        thisText = templateIFCA
        thisText = thisText.replace('NAMEOFJOB', jobname)
    
    thisText = thisText.replace('CMSSWRELEASE', cmsswrelease)
    thisText = thisText.replace('CARD', cardLocation)
    thisText = thisText.replace('EXENAME', exeName)
    thisText = thisText.replace('OUTPUT', outputfilename)
    thisText = thisText.replace('FIRSTEVENT', str(initevent))
    thisText = thisText.replace('LASTEVENT', str(endevent))
    thisText = thisText.replace('INPUT', fileName)
    
    routput = open('submit_{0}.sh'.format(jobname), 'w')
    routput.write(thisText)
    routput.close()
    
    script = open(launcher, 'a')
    script.write('chmod +x submit_{0}.sh\n'.format(jobname))
    if logfilename == 'none' and errfilename == 'none':
        if site == 'ifca':
            script.write('qsub submit_{0}.sh\n'.format(jobname))
        if site == 'lsflxplus':
            script.write('bsub -q {0} submit_{1}.sh\n'.format(queue, jobname))
    else:
        if site == 'ifca':
            script.write('qsub -o {0} -e {1} submit_{2}.sh\n'.format(logfilename, errfilename, jobname))
        if site == 'lsflxplus':
            script.write('bsub -o {0} -e {1} -q {2} submit_{3}.sh\n'.format(logfilename, errfilename, queue, jobname))

    script.close()


##############################################################################################################################################
def prepareJobs(exeName, inputPath, fileName, numberOfEvents, outputDirectory, logLocation, eventsPerJob, cardLocation, cmsswrelease, launcher, queue, site):

    nameTemplate = fileName[0:fileName.find('.root')]
    chunkCounter = 0
    for ji in range(0, numberOfEvents, eventsPerJob):
        jobname = '{0}_chunk{1}'.format(nameTemplate, str(chunkCounter))
        inputfilename = '{0}/{1}.root'.format(inputPath, nameTemplate)
        outputfilename = '{0}/{1}_chunk{2}.root'.format(outputDirectory, nameTemplate, str(chunkCounter))
        initevent = chunkCounter * eventsPerJob
        endevent = (chunkCounter + 1) * eventsPerJob if (chunkCounter + 1) * eventsPerJob < numberOfEvents else numberOfEvents
        if logLocation == 'none':
            prepareJob(exeName, inputfilename, initevent, endevent, outputfilename, 'none', 'none', jobname, cardLocation, cmsswrelease, launcher, queue, site)
        else:
            logfilename = '{0}/{1}_chunk{2}.log'.format(logLocation, nameTemplate, str(chunkCounter))
            errfilename = '{0}/{1}_chunk{2}.err'.format(logLocation, nameTemplate, str(chunkCounter))
            prepareJob(exeName, inputfilename, initevent, endevent, outputfilename, logfilename, errfilename, jobname, cardLocation, cmsswrelease, launcher, queue, site)
        chunkCounter = chunkCounter + 1




##########################################################################################
if __name__ == '__main__':

    parser = optparse.OptionParser(usage='usage: %prog [options] path', version='%prog 1.0')
    parser.add_option('-n', '--number'   , action='store', type='string', dest='eventsPerJob',    default='',     help='Number of events per job.')
    parser.add_option('-o', '--output'   , action='store', type='string', dest='outputDirectory', default='', help='Output directory.')
    parser.add_option('-c', '--card'     , action='store', type='string', dest='cardLocation',    default='', help='Location of delphes card.')
    parser.add_option('-l', '--logs'     , action='store', type='string', dest='logLocation',     default='none',     help='Location of the logs.')
    parser.add_option('-r', '--release'  , action='store', type='string', dest='cmsswrelease',    default='',     help='Location of CMSSW_X_Y_Z')
    parser.add_option('-m', '--launcher' , action='store', type='string', dest='launcher',        default='',     help='Name of launcherfile.')
    parser.add_option('-q', '--queue'    , action='store', type='string', dest='queue',           default='8nh',     help='Queue name.')
    parser.add_option('-s', '--site'     , action='store', type='string', dest='site',            default='lsf',     help='Running site: ifca, lsflxplus, condorlxplus.')
    (opts, args) = parser.parse_args()

    if not opts.launcher or not opts.eventsPerJob or not opts.outputDirectory or not opts.cardLocation or not opts.cmsswrelease or len(args) < 1 or not opts.site:
        print('Error: You need to specify --launcher name --number NumberOfEventsPerJob --output OutputDirectory --card CardLocation --release cmsswrelease --site theSite inputPath')
        sys.exit()

    logLocation = opts.logLocation
    queue = opts.queue
    site = opts.site

    launcher = opts.launcher
    cardLocation = opts.cardLocation
    eventsPerJob = int(opts.eventsPerJob)
    outputDirectory = opts.outputDirectory
    cmsswrelease = opts.cmsswrelease 
    inputPath = args[0]
    
    if site != 'ifca' and site != 'lsflxplus' and site != 'condorlxplus':
        print('Error: The site is not valid.')
        sys.exit()

    if not os.path.isdir(inputPath):
        print('Error: The input path was not found.')
        sys.exit()
    
    if os.path.isfile(launcher):
        print('Error: Launcher already exists.')
        sys.exit()
  
    exeName = '{}/DelphesCMSFWLite'.format(os.getcwd())
    if not os.path.isfile(exeName):
        print('Error: No DelphesCMSSFWLite detected in the working directory.')
        sys.exit()

    listOfFiles = []
    for ji in os.listdir(inputPath):
        if ji.find('.root') == -1:
            continue
        listOfFiles.append(ji)
    if len(listOfFiles) == 0:
        print('Error: No root files were found in the input path')
        sys.exit()
    
    if not os.path.isfile(cardLocation):
        print('Error: Delphes card was not found.')
        sys.exit()

    if eventsPerJob < 1:
        print('Error: Number of events per job should be greater than 0')
        sys.exit()
    
    if not os.path.isdir(outputDirectory):
        os.mkdir(outputDirectory)
    
    if logLocation != 'none' and not os.path.isdir(logLocation):
        os.mkdir(logLocation)

    
    jcount = 0
    for ji in listOfFiles:
        thefile = r.TFile.Open('{0}/{1}'.format(inputPath,ji))
        tree = thefile.Get("Events")
        if not tree:
            print('Error: File {0} is not a valid GEN-SIM file.'.format(ji))
            continue
        numberOfEvents = tree.GetEntries()
        thefile.Close()
        prepareJobs(exeName, inputPath, ji, numberOfEvents, outputDirectory, logLocation, eventsPerJob, cardLocation, cmsswrelease, launcher, queue, site)
        jcount = jcount + 1 

  



