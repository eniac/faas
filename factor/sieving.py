import math
from math import log
import os
from os import path
import time
from time import sleep
import shutil
import sys
import subprocess
import gzip
import re
from collections import deque
import threading

import utils
import filtering
import cadoprograms
import cadoparams
import cadotask
from cadotask import Polynomials

import logging
logger = logging.getLogger('Sieve')
cmd_logger = logging.getLogger('cmd')


class Sieve:

    def __init__(self, parameters, poly_file):
        self.parameters = parameters
        params = self.parameters.myparams({
            "name": str,
            "workdir": str,
            "batch_size": int
        }, ['tasks', 'sieve'])
        self.name = params.get("name")
        self.workdir = params.get("workdir")
        self.batch_size = params.get("batch_size")

        # variables that are accessed from multiple threads
        self.stage = 'sieve' 
        self.stage_lock = threading.Lock()
        self.finished = False
        self.finished_lock = threading.Lock()
        self.queue = deque()
        self.queue_lock = threading.Lock()
        self.rels_total = 0
        self.rels_total_lock = threading.Lock()
        # if rels_wanted is not set, then specify a default initial value based on the large prime bounds
        self.rels_wanted = parameters.myparams({"rels_wanted": 0}, ["tasks", "sieve", "sieving", "las"]).get("rels_wanted")
        if self.rels_wanted == 0:
            # taking into account duplicates, the initial value
            # pi(2^lpbr) + pi(2^lpba) should be good
            paths = ["tasks", "sieve", "sieving", "las"]
            nr = 2 ** self.parameters.myparams({"lpbr": int}, paths).get("lpbr")
            na = 2 ** self.parameters.myparams({"lpba": int}, paths).get("lpba")
            nra = int(nr / log(nr) + na / log(na))
            self.rels_wanted = nra
        self.rels_wanted_lock = threading.Lock()
        
        # compile regular expressions here for speed
        self.relation_re = re.compile("(-?\d*),(\d*):(.*)")
        self.relation_total_re = re.compile("# Total (\d+) reports")
        self.relation_file_re = re.compile("%s[.]sieving[.](\d+)-(\d+)[.]gz" % self.name)

        self.completed_factorization = False

        self.poly_file = poly_file
        self.poly = None
        with open(self.poly_file, 'r') as f:
            self.poly = Polynomials(f.readlines())
        self.reldir = os.path.join(self.workdir, self.name + ".upload")
        self.msieve_dat_file = os.path.join(self.workdir, self.name + ".dat")

        # The relation files that we will pass on to filtering
        # TODO: we might not need this variable
        self.relation_files = []
        # The set of files that we have already seen and processed
        self.seen_files = set()

        self.fb_paths = None
        self.freerel_output = None
        self.generator = None
        self.start_time = None

    def generate_sieving_task_commands(self):
        paths = ["tasks", "sieve", "sieving", "las"]
        program_to_run = cadoprograms.Las
            
        progparams = self.parameters.myparams(program_to_run.get_accepted_keys(), paths)
        progparams.pop("q0", None)
        progparams.pop("q1", None)
        progparams.pop("factorbase", None)
        progparams.pop("out", None)
        progparams.pop("poly", None)
        progparams.pop("stats-stderr", None)

        params = self.parameters.myparams({
            "alim": int,
            "rlim": int,
            "qrange": int,
            "lpba": int,
            "qmin": [int],
            "execpath": str,
        }, paths)

        start = 0
        if "qmin" in params:
            start = params.get("qmin")
        else:
            start = params.get("alim")

        last_q1 = start
        qrange = params.get("qrange")
        qmax = 2 ** params.get("lpba")

        while (True):
            q0 = last_q1
            q1 = q0 + qrange
            last_q1 = q1
            # the special q values should always remain below 2^lpba
            if (q0 >= qmax):
                logger.debug("Stopping attempt to create tasks with special q values greater than 2**lpba")
                break
            # allow q1 to be at most the algebraic large prime bound
            q1 = min(q1, qmax)

            out_file = os.path.join(self.reldir, self.name + '.sieving.' + str(q0) + '-' + str(q1) + '.gz')

            # check if out_file is one of our previously seen files
            if out_file in self.seen_files:
                continue

            if len(self.fb_paths) > 1:  # is twoalgosides
                program = cadoprograms.Las(q0=q0, q1=q1, factorbase0=self.fb_paths[0], factorbase1=self.fb_paths[1], out=out_file, poly=self.poly_file, stats_stderr=True, **progparams)
            else:
                program = cadoprograms.Las(q0=q0, q1=q1, factorbase=self.fb_paths[0], out=out_file, poly=self.poly_file, stats_stderr=True, **progparams)
            cmd = program.make_command_line()
            cmd_logger.debug(cmd)
            yield (out_file, cmd)

    def submit_batch(self):
        jobs = []
        for i in range(self.batch_size):
            try:
                jobs.append(next(self.generator))
            except StopIteration:
                logger.debug("Unable to generate any more tasks")
                # we cannot generate any more commands
                break

        logger.info("Submitting %d additional jobs", len(jobs))

        i = 0
        num_submitted = 0
        batch_file = str(os.path.join(self.workdir, 'sieving.sh'))
        for filen, cmd in jobs:
            with open(batch_file, 'w') as f:
                f.write("#!/bin/sh\n")
                f.write("#SBATCH -p factor\n")
                f.write("#SBATCH -J %s\n" % filen.split('/')[-1])
                f.write("#SBATCH -n 1\n")
                f.write("#SBATCH -c 2\n")
                f.write("#SBATCH -s\n")
                f.write("#SBATCH --requeue\n")
                f.write("#SBATCH --output=/dev/null\n")
                f.write("srun --output=%s.out %s 2>&1\n" % (filen, cmd))
                f.write("wait\n")
            os.chmod(batch_file, 0o755)
            utils.run_command('sbatch ' + batch_file)

            if i >= 100:
                logger.debug("Submitted %d/%d sieve jobs", num_submitted, len(jobs))
                i = 0
                if self.is_finished():
                    return

            num_submitted += 1
            i += 1

    def verify_relation(self, line):
        """ Check that the primes listed for a relation divide the value of
            the polynomials """
        match = self.relation_re.match(line)
        if match:
            a, b, rest = match.groups()
            a, b = int(a), int(b)
            sides = rest.split(":")
            assert len(sides) == 2
            for side, primes_as_str in enumerate(sides):
                value = self.poly.get_polynomial(side).eval_h(a, b)
                primes = [int(s, 16) for s in primes_as_str.split(",")]
                for prime in primes:
                    if value % prime != 0:
                        return False
            return True
        return None

    def process_relation_file(self, filename):
        if not os.path.isfile(filename):
            logger.warning("File '%s' does not exist", filename)
            return 0

        try:
            relations = []
            count = 0
            # check some relations in the file before importing
            relations_to_check = 10
            with gzip.open(filename, 'rt', encoding='utf-8') as f:
                for line in f:
                    if relations_to_check > 0:
                        result = self.verify_relation(line)
                        if result is True:
                            relations_to_check -= 1
                        elif result is False:
                            return 0
                        else:  # Did not match: try again
                            pass

                    if line[0] != '#':
                        relations.append(line)
                        count += 1
            if relations_to_check == 0:
                with open(self.msieve_dat_file, 'at', encoding='utf-8') as f:
                    f.write(''.join(relations))

            # return the relation count in this file
            return count

        except Exception as e:
            logger.warning("Exception in file '%s': %s", filename, e)

        # This will be reached on an exception or if at least 10 relations do not check out
        return 0

    def check_relation_file(self, filename):
        # check that the .out file exists (otherwise, the .gz file might not be complete)
        outfile = filename + '.out'
        outfile_count = 0
        relfile_count = 0
        match = None
        if not os.path.isfile(outfile):
            return 0

        outfile_count = 0
        with open(outfile, 'rt', encoding='utf-8') as f:
            # Check for line matching '# Total 12377 reports [0.00429s/r, 18.0r/sq]'
            match = self.relation_total_re.findall(f.read())
            if match:
                outfile_count = int(match[0])

        return outfile_count

    def import_relations(self, location):
        gz_files = []
        if os.path.isfile(location):
            with open(location, 'r') as f:
                gz_files = [i for i in f if self.relation_file_re.match(i) ]
        elif os.path.isdir(location):
            d = os.listdir(location)
            all_gz_files = [ str(os.path.join(location, i)) for i in d if self.relation_file_re.match(i) ]
            unseen_gz_files = [ i for i in all_gz_files if i not in self.seen_files ]
            gz_files = [ i for i in unseen_gz_files if self.check_relation_file(i) ]
        self.seen_files.update(set(gz_files))
        return gz_files

    def set_stage(self, stage):
        with self.stage_lock:
            self.stage = stage

    def get_stage(self):
        with self.stage_lock:
            return self.stage

    def set_finished(self):
        with self.finished_lock:
            self.finished = True

    def is_finished(self):
        with self.finished_lock:
            return self.finished
    
    def run_slurm_thread(self):
        while not self.is_finished():
            # check if there are enough jobs running
            jobs_out = utils.run_command("squeue -t PENDING,RUNNING,COMPLETING").strip().split('\n')
            logger.info("Number of queued jobs: %d", len(jobs_out) - 1)
            if len(jobs_out) < 2 * self.batch_size:
                self.submit_batch()
            sleep(10)

    def run_sieving(self):    
        logger.info("Starting Sieving...")
        self.start_time = time.time()

        if os.path.isfile(self.msieve_dat_file):
            logger.info("Removing existing .dat file for msieve")
            os.remove(self.msieve_dat_file)

            with open(self.msieve_dat_file, 'wt', encoding='utf-8') as f:
                f.write(str(self.parameters.myparams({'N': int}, [])['N']))
                f.write('\n')

        # check if we should import relations
        import_relation_file = self.parameters.myparams({'import': None}, ['tasks', 'sieve']).get('import')
        if import_relation_file:
            logger.info("Importing relations from file '%s'", import_relation_file)
            imported_files = self.import_relations(import_relation_file)
            self.queue_extend(imported_files)
            logger.info("Found %d relation files in file '%s'", len(imported_files), self.reldir) 

        if not os.path.exists(self.reldir):
            # create the directory for relations if it does not yet exist
            logger.info("Creating directory for relations %s", self.reldir)
            os.makedirs(self.reldir)
        else:
            # check if there are relations in the directory already and add them to the import queue.
            # We do this check once outside the loop so that we know which tasks not to regenerate
            logger.info("Importing relations files from directory '%s;", self.reldir)
            imported_files = self.import_relations(self.reldir)
            self.queue_extend(imported_files)
            logger.info("Found %d relation files in directory '%s'", len(imported_files), self.reldir) 

        # generate sieving task commands, skipping files that have already been generated
        self.generator = self.generate_sieving_task_commands()

        # spawn a thread to launch sieving tasks until sieving is finished
        slurm_thread = threading.Thread(target=self.run_slurm_thread)
        slurm_thread.start()

        # spawn a thread to print out a status message periodically
        status_thread = threading.Thread(target=self.run_status_thread)
        status_thread.start()

        # spawn a thread to filter periodically until sieving is finished
        filter_thread = threading.Thread(target=self.run_filter_thread)
        filter_thread.start()

        # wait until filter thread completes
        try:
            filter_thread.join()
            slurm_thread.join()
            status_thread.join()
        finally:
            utils.run_command("scancel -p factor")

    def queue_empty(self):
        with self.queue_lock:
            return False if self.queue else True 

    def queue_pop(self):
        with self.queue_lock:
            return self.queue.popleft()

    def queue_extend(self, elements):
        with self.queue_lock:
            self.queue.extend(elements)

    def get_rels_total(self):
        with self.rels_total_lock:
            return self.rels_total 

    def set_rels_total(self, value):
        with self.rels_total_lock:
            self.rels_total = value

    def get_rels_wanted(self):
        with self.rels_wanted_lock:
            return self.rels_wanted 

    def set_rels_wanted(self, value):
        with self.rels_wanted_lock:
            self.rels_wanted = value

    def run_filter_thread(self):
        # first, initialize filtering state
        filtering.init(self.parameters, self.freerel_output)

        # a list of new files to be passed in to filtering
        new_files = []
        while True:
            # check if there are relations in the directory already and add them to the import deque
            imported_files = self.import_relations(self.reldir)
            self.queue_extend(imported_files)

            # process files from queue until we reach rels_wanted or the queue is empty
            while not self.queue_empty():
                filename = self.queue_pop()
                count = self.process_relation_file(filename)
                self.set_rels_total(self.get_rels_total() + count)

                self.relation_files.append(filename)
                new_files.append(filename)
                logger.info("Found %d new relations in %s. Total %d/%d (%.2f%%)",\
                            count, filename, self.get_rels_total(), self.get_rels_wanted(), 100 * self.get_rels_total() / self.get_rels_wanted())
                if self.get_rels_total() >= self.get_rels_wanted():
                    logger.info("Reached rels_wanted with %d/%d relations", self.get_rels_total(), self.get_rels_wanted())
                    # stop processing output files when we've reached our target
                    break

            if self.get_rels_total() >= self.get_rels_wanted():
                self.set_stage('filter')
                rels_additional = filtering.run(new_files)
                new_files = []
                if rels_additional == 0:
                    break
                elif rels_additional == -1:
                    rels_addtional = math.ceil(self.get_rels_wanted() * 0.1)
                elif rels_additional == -2:
                    self.completed_factorization = True
                    break
                elif rels_additional > 0:
                    self.set_rels_wanted(self.get_rels_wanted() +  rels_additional)
                else:
                    logger.error('unexpected value for rels_additional: ', rels_additional)
                    sys.exit(1)
            else:
                self.set_stage('sieve')
            sleep(10)
        self.set_finished()

    def run_status_thread(self):
        while not self.is_finished():

            # give the user a status update
            stage = self.get_stage()
            if stage == 'sieve':
                elapsed = time.time() - self.start_time
                rels_total = self.get_rels_total()
                rels_wanted = self.get_rels_wanted()

                rate = rels_total / elapsed
                eta = 0
                if rate > 0:
                    eta = (rels_wanted - rels_total) / rate
                    logger.info("Status: %d/%d relations at %d rels/sec - elapsed: %s, ETA: %s", rels_total, rels_wanted, int(rate), utils.str_time(elapsed), utils.str_time(eta))
            elif stage == 'filter':
                logger.info("Status: performing filtering")
                    
            sleep(10)
    
    def run_factor_base(self):
        logger.info("Starting FactorBase...")
        paths = ["tasks", "sieve", "factorbase", "makefb"]
        fb_params = self.parameters.myparams({"gzip": True, "I": int, "rlim": int, "alim": int}, paths)

        program_to_run = cadoprograms.MakeFB
        progparams = self.parameters.myparams(program_to_run.get_accepted_keys(), paths)
        progparams.setdefault("maxbits", fb_params["I"] - 1)
        input_files = {"poly": self.poly_file}
        merged_args = dict(progparams.items() | input_files.items())

        use_gz = ".gz" if fb_params["gzip"] else ""

        twoalgsides = self.poly.polyg.degree > 1
        outpath = os.path.join(self.workdir, "factorbase")
        if not os.path.exists(outpath):
            logger.info("Creating directory for factorbase %s", outpath)
            os.makedirs(outpath)
        outputfiles = ["%s/roots%d%s" % (outpath, i, use_gz) for i in range(1 + twoalgsides)]
        if not twoalgsides:
            p = cadoprograms.MakeFB(out=outputfiles[0], lim=fb_params["alim"], **merged_args)
            command_line = p.make_command_line()
            cmd_logger.debug(command_line)
            utils.run_command(command_line)
            logger.info("FactorBase: Created outputfile %s", outputfiles[0])
        else:
            p = cadoprograms.MakeFB(out=outputfiles[1], side=0, lim=fb_params["rlim"], **merged_args)
            command_line = p.make_command_line()
            cmd_logger.debug(command_line)
            utils.run_command(command_line)
            logger.info("FactorBase: Created outputfile %s", outputfiles[1])

            p = cadoprograms.MakeFB(out=outputfiles[0], side=1, lim=self.params["alim"], **merged_args)
            command_line = p.make_command_line()
            cmd_logger.debug(command_line)
            utils.run_command(command_line)
            logger.info("FactorBase: Created outputfile %s", outputfiles[0])
        self.fb_paths = outputfiles

    def parse_freerel_output(self, text):
        wanted_regex = {
            'nfree': (r'# Free relations: (\d+)', int),
            'nprimes': (r'Renumbering struct: nprimes=(\d+)', int)
        }
        found = {}
        for line in text.split("\n"):
            for (key, (regex, datatype)) in wanted_regex.items():
                match = re.match(regex, line)
                if match:
                    found[key] = datatype(match.group(1))
        return found

    def run_freerel(self):
        logger.info("Starting FreeRel...")
        paths = ["tasks", "sieve", "freerel"]

        program_to_run = cadoprograms.FreeRel
        progparams = self.parameters.myparams(program_to_run.get_accepted_keys(), paths)
        progparams.setdefault("addfullcol", True)
        input_files = {"poly": self.poly_file}
        merged_args = dict(progparams.items() | input_files.items())

        use_gz = ".gz" if self.parameters.myparams({"gzip": True}, paths) else ""

        outpath = os.path.join(self.workdir, "freerel")
        if not os.path.exists(outpath):
            logger.info("Creating directory for freerel %s", outpath)
            os.makedirs(outpath)
        freerelfilename = "%s/freerel%s" % (outpath, use_gz)
        renumberfilename = "%s/renumber%s" % (outpath, use_gz)
        p = cadoprograms.FreeRel(renumber=renumberfilename, out=freerelfilename, **merged_args)
        command_line = p.make_command_line()
        cmd_logger.debug(command_line)
        stdout, stderr = utils.run_command(command_line, include_stdout=True, include_stderr=True, logger=logger)
        freerel_output = self.parse_freerel_output(stderr)
        freerel_output["freerelfilename"] = freerelfilename
        freerel_output["renumberfilename"] = renumberfilename
        logger.info("FreeRel: Created outputfile %s", freerelfilename)
        logger.info("FreeRel: Created outputfile %s", renumberfilename)
        self.freerel_output = freerel_output

    def run(self):
        self.run_factor_base()
        self.run_freerel()
        self.run_sieving()
        return self.relation_files

