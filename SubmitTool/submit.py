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

template = """
#!/bin/bash
pushd CMSSWRELEASE/src
eval `scramv1 runtime -sh`
pushd
./DelphesCMSFWLite CARD OUTPUT FIRSTEVENT LASTEVENT INPUT 
"""

##############################################################################################################################################
def prepareJob(fileName, initevent, endevent, outputfilename, logfilename, errfilename, jobname, cardLocation, cmsswrelease, launcher, queue):

    thisText = template
    thisText = thisText.replace('CMSSWRELEASE', cmsswrelease)
    thisText = thisText.replace('CARD', cardLocation)
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
        script.write('qsub -q {0} submit_{1}.sh\n'.format(queue, jobname))
    else:
        script.write('qsub -o {0} -e {1} -q {2} submit_{3}.sh\n'.format(logfilename, errfilename, queue, jobname))
    script.close()


##############################################################################################################################################
def prepareJobs(inputPath, fileName, numberOfEvents, outputDirectory, logLocation, eventsPerJob, cardLocation, cmsswrelease, launcher, queue):

    nameTemplate = fileName[0:fileName.find('.root')]
    chunkCounter = 0
    for ji in range(0, numberOfEvents, eventsPerJob):
        jobname = '{0}_chunk{1}'.format(nameTemplate, str(chunkCounter))
        inputfilename = '{0}/{1}.root'.format(inputPath, nameTemplate)
        outputfilename = '{0}/{1}_chunk{2}.root'.format(outputDirectory, nameTemplate, str(chunkCounter))
        initevent = chunkCounter * eventsPerJob
        endevent = (chunkCounter + 1) * eventsPerJob if (chunkCounter + 1) * eventsPerJob < numberOfEvents else numberOfEvents
        if logLocation == 'none':
            prepareJob(inputfilename, initevent, endevent, outputfilename, 'none', 'none', jobname, cardLocation, cmsswrelease, launcher, queue)
        else:
            logfilename = '{0}/{1}_chunk{2}.log'.format(logLocation, nameTemplate, str(chunkCounter))
            errfilename = '{0}/{1}_chunk{2}.err'.format(logLocation, nameTemplate, str(chunkCounter))
            prepareJob(inputfilename, initevent, endevent, outputfilename, logfilename, errfilename, jobname, cardLocation, cmsswrelease, launcher, queue)
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
    (opts, args) = parser.parse_args()

    if not opts.launcher or not opts.eventsPerJob or not opts.outputDirectory or not opts.cardLocation or not opts.cmsswrelease or len(args) < 1:
        print('Error: You need to specify --launcher name --number NumberOfEventsPerJob --output OutputDirectory --card CardLocation --release cmsswrelease inputPath')
        sys.exit()

    logLocation = opts.logLocation
    queue = opts.queue

    launcher = opts.launcher
    cardLocation = opts.cardLocation
    eventsPerJob = int(opts.eventsPerJob)
    outputDirectory = opts.outputDirectory
    cmsswrelease = opts.cmsswrelease 
    inputPath = args[0]
    
    if not os.path.isdir(inputPath):
        print('Error: The input path was not found.')
        sys.exit()
    
    if os.path.isfile(launcher):
        print('Error: Launcher already exists.')
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
        prepareJobs(inputPath, ji, numberOfEvents, outputDirectory, logLocation, eventsPerJob, cardLocation, cmsswrelease, launcher, queue)
        jcount = jcount + 1 

  



