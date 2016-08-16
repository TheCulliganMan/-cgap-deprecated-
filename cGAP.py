#!/usr/bin/env python
import argparse
import sys
import os
import subprocess as sp
import time
import threading
import shutil
import shlex
from collections import deque
from getBlastHits import get_blast_hits
try:
    import Queue
except ImportError:
    import queue as Queue


"""
                     ,ad8888ba,        db        88888888ba
                    d8"'    `"8b      d88b       88      "8b
                   d8'               d8'`8b      88      ,8P
         ,adPPYba, 88               d8'  `8b     88aaaaaa8P'
        a8"     "" 88      88888   d8YaaaaY8b    88""""""'
        8b         Y8,        88  d8""""""""8b   88
        "8a,   ,aa  Y8a.    .a88 d8'        `8b  88
         `"Ybbd8"'   `"Y88888P" d8'          `8b 88
        Consensus Genomic Alignment Program

#Ryan Culligan
#Julian Egger
#Richard Bankoff
"""


class CGAP():
    # Main function used to take in arguments and check for errors.
    def __init__(self):
        """
        Inputs:
            SELF
        Returns:
            None
        Description:
            Parses arguments.
            Initializes variables.
            Creates commandsLists.
            A little too much for most classes, but helps keep the variables under control.
        """

        parser = argparse.ArgumentParser(description='Start cGAP: Script for running Consensus Gene Assembly Program')

        parser.add_argument('-refs_path', action="store", dest="directory", help="The path where cGAP gene references are stored.")  # Directory of query files
        parser.add_argument('-forward', action="store", dest="forwardReads", nargs="+", help="The FORWARD fastq samples that will be used for the run.")  # Forward Short Read Files
        parser.add_argument('-reverse', action="store", dest="reverseReads", nargs="+", help="The REVERSE fastq samples that will be used for this run.")  # Reverse Short Read Files
        parser.add_argument('-samples', action="store", dest="sampleNames", nargs="+", help="The the sample name you would like sequences to be named after.")  # Sample Names
        parser.add_argument('-work_path', action="store", dest="workPath", help="The full path directory in which cGAP files are located and will be created.")  # Full path for working directory
        parser.add_argument('-formatdb', action="store_true", dest="formatdb", help="Turns on SRDB format for cGAP samples")  # Option for formatting database or not
        parser.add_argument('-c', action="store", dest="cores", type=int, help="The number of tasks cGAP should run at a time.  Defaults to 1.")  # gets number of cores
        parser.add_argument('-maskDepth', action="store", dest="maskDepth", type=int, help="The depth to which sequences should be masked.  Defaults to 4.")

        exCommand = "Example: python cGAP.py -refs_path references -forward mery_3a1.fastq -reverse mery_3a2.fastq -samples mery3 -work_path /cGAP -formatdb -c 32"

        if len(sys.argv) < 2:
            parser.print_usage()
            print(exCommand)
            sys.exit(1)

        args = parser.parse_args()

        error = False
        # Check if all arguments have been entered
        if args.directory is None:
            error = True
            print("\nNo directory for query files entered!!!\n")
        if args.forwardReads is None:
            error = True
            print("\nNo forward reads files entered!!!\n")
        if args.reverseReads is None:
            error = True
            print("\nNo reverse reads files enetered!!!\n")
        if args.sampleNames is None:
            error = True
            print("\nNo sample names entered!!!\n")
        if args.workPath is None:
            error = True
            print("\nNo path to working directory entered!!!\n")

        if args.cores is None:  # default to one core
            args.cores = 1
        if args.maskDepth is None:  # default to mask depth 4
            args.maskDepth = 4

        if error:
            parser.print_usage()
            print(exCommand)
            sys.exit(1)

        self.cores = args.cores
        self.maskDepth = args.maskDepth
        self.directory = args.directory
        self.forwardReads = args.forwardReads
        self.reverseReads = args.reverseReads
        self.sampleNames = args.sampleNames
        self.workPath = args.workPath
        self.formatDB = args.formatdb

        self.geneFileList = []  # list of query files
        self.geneNameList = []  # list of query names taken from file names

        """Sets to ensure no duplicated commands"""
        self.cnsCommands = set([])
        self.maskCommands = set([])

        self.indexCommandList = []
        self.pipelineCommands = []

        # get names of query genes from query file names
        # name is determined as the group of characters before first "." in file name
        path = str(self.directory)
        for filename in os.listdir(path):
            if filename.endswith(".fa") or filename.endswith(".fasta"):
                self.geneFileList.append("{self.directory}/{filename}".format(**locals()))
                self.geneNameList.append(filename.split(".")[0])

        self.commandDict = {geneName: {s: [] for s in args.sampleNames} for geneName in self.geneNameList}
        self.hit_command_dict = {}
        #print self.c

    def getFormatDB(self):
        """
        Inputs:
            SELF
        Returns:
            True
        Description:
            Makes commands to create a blast database.
        """

        allReads = self.forwardReads + self.reverseReads
        formatCommands = []

        for fileName in allReads:
            print ('Formatting srdb for ', fileName)
            formatCommands.append("perl bin/formatSRDB.pl {fileName} ".format(**locals()))

        return formatCommands

    # Function used to create directories for storing and organizing all files created by cGAP
    def makeSampleDirectories(self):
        print ("RUN SAMPLE DIRS")
        """
        Inputs:
            SELF
        Returns:
            True
        Description:
            Makes folders for cGAP to run.
        """

        for s in self.sampleNames:  # make sample directories
            directories = ["{self.workPath}/{s}".format(**locals()),
                           "{self.workPath}/{s}/cnsSeqs.{s}".format(**locals()),
                           "{self.workPath}/{s}/bamFiles.{s}".format(**locals()),
                           "{self.workPath}/fr_readHits".format(**locals()),
                           "{self.workPath}/queryWork".format(**locals())]
            for directory in directories:
                if not os.path.exists(directory):
                    os.makedirs(directory)

        for s, f, r in zip(self.sampleNames, self.forwardReads, self.reverseReads):
            # Create subdirectory in sample directory named after forwardshort reads file w/o extensions
            for directory in (f.split(".fastq")[0], r.split(".fastq")[0]):
                directory = "{self.workPath}/{s}/{directory}".format(**locals())
                if not os.path.exists(directory):
                    os.makedirs(directory)

        return True

    def makeBlastCommands(self, s, q, n, read, searchFastq):
        """
        Inputs:
            SELF,
            s: The sample name that the command is working on.
            n: The name of the gene that the command is working on for the aforementioned sample.
            read: FORWARD exclusive or REVERSE read folder
            searchFastq: the forward or reverse blastdb file to search
        Returns:
            True
        Description:
            Creates commands for each sample / gene / direction to form the .blast files.
        """

        blastFile = '{self.workPath}/{s}/{read}/{read}.{n}.blast'.format(**locals())
        blastDict = q

        referenceFile = "{searchFastq}".format(**locals())

        cmd = 'blastall -p blastn -d {referenceFile} -i {blastDict} -v 500000 -m 8 > {blastFile};   '.format(**locals())
        self.commandDict[n][s] += [cmd]
        return True

    def makeGetHits(self, s, n, read):
        """
        Inputs:
            SELF,
            s: The sample name that the command is working on.
            n: The name of the gene that the command is working on for the aforementioned sample.
            read: FORWARD exclusive or REVERSE
        Returns:
            True
        Description:
            Creates commands for each sample / gene / direction to form the .hits files.
        """

        blastFile = '{self.workPath}/{s}/{read}/{read}.{n}.blast'.format(**locals())  # output file from blast
        cmd = "python bin/getHits.py -blastFile {blastFile} -n 50; ".format(**locals())
        self.commandDict[n][s] += [cmd]
        return True

    def makeHitsCommands(self, s, n, read, fastq_name, hitsDirections):
        """
        Inputs:
            SELF,
            s: The sample name that the command is working on.
            n: The name of the gene that the command is working on for the aforementioned sample.
            read: FORWARD exclusive or REVERSE
            hitsDirections: blasthits exclusive or hitpairs
        Returns:
            True
        Description:
            Creates commands for each sample / gene / direction to form the blasthits/hitpars fastq.
        """
        #if the gene name isn't in the hit command dict add it
        if fastq_name not in self.hit_command_dict:
            self.hit_command_dict[fastq_name] = {"sample": [], "hit_files": [], "gene": []}
        for direction in hitsDirections:
            direction = direction.rsplit(".fastq")[0]
            self.hit_command_dict[fastq_name]['sample'].append('{s}'.format(**locals()))
            self.hit_command_dict[fastq_name]['gene'].append('{n}'.format(**locals()))
            self.hit_command_dict[fastq_name]['hit_files'].append('{self.workPath}/{s}/{direction}/{direction}.{n}.hits'.format(**locals()))

        return True

    def make_hit_list_dictionary(self):
        commands = []
        for output_file, fastq_file in enumerate(self.hit_command_dict):
            hit_files_list = self.hit_command_dict[fastq_file]
            data = {}
            data["genes"] = hit_files_list["gene"]
            data["samples"] = hit_files_list["sample"]
            data["hit_files"] = hit_files_list['hit_files']
            data["fastq_file"] = "{self.workPath}/{fastq_file}".format(**locals())
            data["work_path"] = "{self.workPath}".format(**locals())
            data["direction"] = "{}".format(output_file)
            commands.append(data)
        return commands

    def make_hit_command_items(self):
        """
        Inputs:
            SELF,
            s: The sample name that the command is working on.
            n: The name of the gene that the command is working on for the aforementioned sample.
            read: FORWARD exclusive or REVERSE specification.
        Returns:
            True
        Description:
            Creates commands for each sample / gene reference to for the consensus sequence.
            appends commands to the self.commandDict
        """
        commands = []
        for output_file, fastq_file in enumerate(self.hit_command_dict):
            sample = fastq_file.rsplit(".")[0]

            hit_files_list = self.hit_command_dict[fastq_file]
            genes = "".join((i+" " for i in hit_files_list["gene"]))
            samples = "".join((i+" " for i in hit_files_list["sample"]))
            hit_files_str = "".join((i+" " for i in hit_files_list['hit_files']))

            cmd = "pypy bin/getBlastHits.py "  # THIS WILL NEED FIXING
            cmd += "-fastq_file {self.workPath}/{fastq_file} ".format(**locals())
            cmd += "-work_path {self.workPath} ".format(**locals())
            cmd += "-direction {} ".format(output_file)
            cmd += "-genes {} ".format(genes)
            cmd += "-samples {} ".format(samples)
            cmd += "-hit_files {hit_files_str} ".format(**locals())
            commands += [cmd]
        return commands

    def createQueryFiles(self, q, n):
        """
        Inputs:
            SELF,
            q: The reference file name.
            n: The name of the gene that the command is working on for the aforementioned sample.
        Returns:
            True
        Description:
            Creates commands to copy reference files to a query folder and indexes the fasta
            files for samtools.
        """

        file_name = "{self.workPath}/queryWork/{n}.fa".format(**locals())

        if not os.path.isfile(file_name):
            shutil.copy(q, file_name)

        cmd = "bwa index {file_name};   ".format(**locals())
        cmd += "samtools faidx {file_name};  ".format(**locals())
        self.indexCommandList.append(cmd)

        return True

    def addCmd(self, n, s, cmd):
        """
        Inputs:
            SELF,
            s: The sample name that the command is working on.
            n: The name of the gene that the command is working on for the aforementioned sample.
            cmd: consensus command string for each gene and sample.
        Returns:
            True
        Description:
            Modifies self dictionaries for the makeConsensus function.  Not totally necessary,
            but avoids a lot of repitition.
        """

        if cmd not in self.cnsCommands:
            self.commandDict[n][s].append(cmd)
            self.cnsCommands.add(cmd)
        return True

    def makeConsensus(self, s, n):
        """
        Inputs:
            SELF,
            s: The sample name that the command is working on.
            n: The name of the gene that the command is working on for the aforementioned sample.
        Returns:
            True
        Description:
            Creates commands for each sample / gene reference to for the consensus sequence.
        """

        refFile = "{self.workPath}/queryWork/{n}.fa".format(**locals())
        all_reads = "{self.workPath}/fr_readHits/{n}.{s}.**.fastq".format(**locals())
        masked = "{self.workPath}/{s}/cnsSeqs.{s}/{n}.masked.fasta".format(**locals())
        bamFileOut = "{self.workPath}/{s}/bamFiles.{s}/{n}.final.bam".format(**locals())  # final bam file
        vcfFileOut = "{self.workPath}/{s}/bamFiles.{s}/{n}.vcf"
        depthFileOut = "{self.workPath}/{s}/bamFiles.{s}/{n}.tdv"

        cmd =  'bwa mem {refFile} {all_reads} | '
        cmd += 'samtools view -Su - | '
        cmd += 'novosort -rd -t . - > {bamFileOut}; '
        self.addCmd(n, s, cmd)
        # add command

        cmd =  'cat {bamFileOut} | '
        cmd += 'samtools mpileup -A -ug -f {refFile} - | '
        cmd += 'bcftools call -c | '
        cmd += 'bgzip > {vcfFileOut};'
        self.addCmd(n, s, cmd)

        cmd = 'tabix {vcfFileOut};'
        self.addCmd(n, s, cmd)

        cmd =  "bcftools filter -i'(%QUAL<20) || (%QUAL==999) || (DP <= 3)' {vcfFileOut} | "
        cmd += "bcftools query -f'%CHROM\t%POS\n > {depthFile}; "
        self.addCmd(n, s, cmd)

        cmd =  "bcftools consensus {vcfFileOut} -f {refFile} -m {depthFile} > {cns_file}; "
        self.addCmd(n, s, cmd)

        return True


    def create_pipes(self, pipeline_step):
        """
        Inputs:
            SELF
        Returns:
            True
        Description:
            This function runs functions to modify class objects in the correct to form command lists
            in sample / gene pair dictionaries that can be run in parallel by the main function.

        """
        if pipeline_step == 0:
            print("BLASTING REFERENCES AGAINST FASTQ FILES")
            for q, n in zip(self.geneFileList, self.geneNameList):
                pass
                self.createQueryFiles(q, n)

            for step in range(2):
                for f, r, s in zip(self.forwardReads, self.reverseReads, self.sampleNames):
                    #get the hits
                    for readDirectory in (f, r):  # forward, reverse for each sample
                        read = readDirectory.split(".fastq")[0]
                        for q, n in zip(self.geneFileList, self.geneNameList):
                            if step == 0:
                                self.makeBlastCommands(s, q, n, read, readDirectory)
                            elif step == 1:
                                self.makeGetHits(s, n, read)

        elif pipeline_step == 1:
            print("GETTING BLAST HITS")
            for f, r, s in zip(self.forwardReads, self.reverseReads, self.sampleNames):
                #get the hits
                for readDirectory in (f, r):  # forward, reverse for each sample
                    read = readDirectory.split(".fastq")[0]
                    for q, n in zip(self.geneFileList, self.geneNameList):
                        if readDirectory == f.split(".fastq")[0]:  # for get hits
                            hitsDirections = (f, r)
                        else:
                            hitsDirections = (r, f)
                        self.makeHitsCommands(s, n, read, readDirectory, hitsDirections)

            return self.make_hit_list_dictionary()

        else:
            """GO PARALLEL AGAIN"""
            print("GENERATING CGAP OUTPUT")
            for f, r, s in zip(self.forwardReads, self.reverseReads, self.sampleNames):
                #get the hits
                for readDirectory in (f, r):  # forward, reverse for each sample
                    read = readDirectory.split(".fastq")[0]
                    for q, n in zip(self.geneFileList, self.geneNameList):
                        self.makeConsensus(s, n)
        return True

    def mergePipes(self):
        """
        Inputs:
            SELF
        Returns:
            True
        Description:
            This function parses self.commandDict to form 2D command lists (self.geneCallList) that
            can be run in parallel by the main function.

        """
        nwcmddict = {}
        for geneName in self.commandDict:
            geneCallList = []
            for sample in self.commandDict[geneName]:
                if sample not in nwcmddict:
                    nwcmddict[sample] = []
                isCns = os.path.isfile("{self.workPath}/{sample}/cnsSeqs.{sample}/{sample}.{geneName}.cns.fasta".format(**locals()))
                if not isCns:
                    geneCommand = self.commandDict[geneName][sample]
                    geneCallList = [command for command in geneCommand]
                    nwcmddict[sample].append(geneCallList)
        new = []
        for sample in nwcmddict:
            print(sample)
            new += nwcmddict[sample]
        self.pipelineCommands = new
        return True


def main():
    """
    Inputs:
        None
    Returns:
        True
    Description:
        runs the CGAP class.  Uses commands from the cGAP class to run in parallel.
    """
    for pipeline_step in range(3):

        runCGAP = CGAP()

        runCGAP.makeSampleDirectories()
        status = runCGAP.create_pipes(pipeline_step=pipeline_step)

        if pipeline_step != 1:
            runCGAP.mergePipes()
        else:
            commandList = status

        indexCommands = runCGAP.indexCommandList
        total_genes = len(indexCommands)
        formatCommands = runCGAP.getFormatDB()

        if pipeline_step != 1:
            commandList = runCGAP.pipelineCommands

        if pipeline_step == 0:
            time.sleep(1)  # make sure everything is actually in the right place
            print("INDEX COMMANDS")
            for item in indexCommands:
                sp.call(item, shell=True)  # NEED TO CLEAN THIS UP
            sys.stdout.flush()

        if runCGAP.formatDB and pipeline_step == 0:
            newCommands = formatCommands[0]  # setup first item
            if len(formatCommands) > 1:
                for command in formatCommands[1:]:  # start with the second item
                    newCommands += "& {command}".format(**locals())
            sp.call("{newCommands};".format(**locals()), shell=True)

        #multicore threading
        if pipeline_step is None:  # so this step involves a ton of IO turn to 1 if you want no parallel
            for sample_num, data in enumerate(commandList):
                genes = data["genes"]
                samples = data["samples"]
                hit_files = data["hit_files"]
                fastq_file = data["fastq_file"]
                work_path = data["work_path"]
                direction = data["direction"]
                print("STARTING", sample_num)
                get_blast_hits(genes, samples, hit_files, fastq_file, work_path, direction)

        else:
            threadList = []

            if pipeline_step == 1:  # so this step involves a ton of IO need to turn down the parallel
                cpu_step_arg = 2  # this can be anything.  I wouldn't recommend more than 3 or 4 for traditional hard drives in ra
            else:
                cpu_step_arg = runCGAP.cores

            for i in range(cpu_step_arg):  # NUMBER CPU
                threadList.append("Thread-{}".format(i+1))

            queueLock = threading.Lock()
            workQueue = Queue.Queue(len(commandList))
            threads = []
            threadID = 1

            exitFlag = 0

            class myThread (threading.Thread):
                def __init__(self, threadID, name, q):
                    threading.Thread.__init__(self)
                    self.threadID = threadID
                    self.name = name
                    self.q = q

                def run(self):
                    process_data(self.name, self.q)
                    print("Exiting " + self.name)

            def process_data(threadName, q):
                while not exitFlag:

                    queueLock.acquire()

                    if not workQueue.empty():
                        data = q.get()
                        queueLock.release()
                        status = True  # NEED TO ADD ERROR CHECKING

                        if pipeline_step != 1:
                            for item in data:
                                status = sp.call(item, shell=True)
                        else:
                            genes = data["genes"]
                            samples = data["samples"]
                            hit_files = data["hit_files"]
                            fastq_file = data["fastq_file"]
                            work_path = data["work_path"]
                            direction = data["direction"]
                            get_blast_hits(genes, samples, hit_files, fastq_file, work_path, direction)

                    else:
                        queueLock.release()

            # Create new threads
            for tName in threadList:
                thread = myThread(threadID, tName, workQueue)
                thread.start()
                threads.append(thread)
                threadID += 1

            # Fill the queue
            queueLock.acquire()
            for command in commandList:
                workQueue.put(command)
            queueLock.release()

            # Wait for queue to empty
            while not workQueue.empty():
                time.sleep(1)

            # Notify threads it's time to exit
            exitFlag = 1

            # Wait for all threads to complete
            for thread in threads:
                thread.join()

            print("Exiting Main Thread")

            exitFlag = 0


if __name__ == "__main__":
    main()
