import os
import time
from time import sleep
import utils
import threading
from cadotask import Polynomials
import cadoprograms
import cadoparams
import sys
import heapq
import subprocess

import logging
logger = logging.getLogger("Polysel")
cmd_logger = logging.getLogger('cmd')

'''
Polynomial Selection
'''

class Polysel:

    def __init__(self, parameters):
        self.parameters = parameters

        params = parameters.myparams({
            "name": str,
            "workdir": str
        }, ["tasks", "polyselect"])
        self.name = params.get("name")
        self.workdir = params.get("workdir")
        self.msieve_fb_path = os.path.join(self.workdir, self.name + ".fb")

    # yield commands one by one
    def generate_polysel1_task_commands(self, params):

        adnext = params["admin"]
        admax = params["admax"]
        adrange = params["adrange"]
        outputdir = params["outputdir"]
        while adnext < admax:
            adstart = adnext
            adend = adstart + adrange
            adend = min(adend, admax) 
            outputfile = str(os.path.join(outputdir, "%s.polyselect1.%d-%d" % (self.name, adstart, adend)))

            # copy and remove illegal Polyselect2l parameters
            progparams = params.copy()
            progparams.pop("admin", None)
            progparams.pop("admax", None)
            progparams.pop("adrange", None)
            progparams.pop("nrkeep", None)
            progparams.pop("batch_size", None)
            progparams.pop("outputdir", None)

            program = cadoprograms.Polyselect2l(admin=adstart, admax=adend, stdout=outputfile, **progparams)

            task_command = program.make_command_line() 
            cmd_logger.debug(task_command)
            adnext = adend
            yield task_command, outputfile

    def submit_polysel1_batch(self, generator, params):
        task_commands = []
        batch_size = params["batch_size"]
        batch_file = str(os.path.join(params["outputdir"], "polysel1.sh"))
        for i in range(batch_size):
            try:
                task_commands.append(next(generator))
            except StopIteration as e:
                logger.debug("Cannot generate any more polysel1 tasks")
                break
        logger.info("Submitting %d additional polysel1 tasks", len(task_commands))

        outputfiles = []
        i = 0
        num_submitted = 0
        for command, outputfile in task_commands:
            with open(batch_file, 'w') as f:
                outputfiles.append(outputfile)
                f.write("#!/bin/sh\n")
                f.write("#SBATCH -p factor\n")
                f.write("#SBATCH -J %s\n" % outputfile.split('/')[-1])
                f.write("#SBATCH -n 1\n")
                f.write("#SBATCH -c 2\n")
                f.write("#SBATCH -s\n")
                f.write("#SBATCH --requeue\n")
                f.write("#SBATCH --output=/dev/null\n")
                f.write("srun -s -c 2 -n 1 -p factor -J polysel1 " + command + " 2>&1 &\n")
                f.write("wait\n")
            os.chmod(batch_file, 0o755)
            # TODO: use the stdout of sbatch to determine job IDs, and only continue when all of the jobs have finished.
            stdout = utils.run_command("sbatch " + batch_file)

            # rate-limit job submission so slurm is not overwhelmed
            if i >= 100:
                logger.debug("Submitted %d/%d polysel1 jobs", num_submitted, len(task_commands))
                #time.sleep(1)
                i = 0
            num_submitted += 1
            i += 1

        while (True):
            jobs_out = utils.run_command("squeue -t PENDING,RUNNING,COMPLETING").strip().split('\n')
            if len(jobs_out) - 1 <= 0:
                break
            logger.info("Number of queued jobs: %d", len(jobs_out) - 1)
            sleep(10)

        polys = []
        for filen in outputfiles:
            with open(filen) as f:
                polys += list(utils.parse_poly(f.read()))
        return polys

    def run_polysel1(self, parameters):
        logger.info("Starting polysel1")
        start_time = time.time()

        # get parameters for polysel1
        polysel1_paths = ["tasks", "polyselect", "polyselect1", "polyselect2l"]
        polysel1_program = cadoprograms.Polyselect2l
        polysel1_keys = polysel1_program.get_accepted_keys()
        polysel1_keys.update({"batch_size": int, "nrkeep": int, "admin": 0, "admax": int, "adrange": int})
        polysel1_params = parameters.myparams(polysel1_keys, polysel1_paths)
        polysel1_params.update({"outputdir": str(os.path.join(self.workdir, self.name + ".upload"))})
        if not os.path.exists(polysel1_params["outputdir"]):
            logger.info("Creating directory for polysel1 output files %s", polysel1_params["outputdir"])
            os.makedirs(polysel1_params["outputdir"])

        nrkeep = polysel1_params["nrkeep"]
        batch_size = polysel1_params["batch_size"]

        generator = self.generate_polysel1_task_commands(polysel1_params)
        polys = []
        while True:
            new_polys = self.submit_polysel1_batch(generator, polysel1_params)
            if len(new_polys) == 0:
                break
            polys += new_polys

        polysel1_bestpolys = heapq.nsmallest(nrkeep, polys, lambda s: s.lognorm)
        if not polysel1_bestpolys:
            logger.critical("No polys recieved from polysel1. Exiting...")
            sys.exit(0)
        
        logger.info("Found best %d polynomials", len(polysel1_bestpolys))
        logger.info("Polysel1 finished in %s", utils.str_time(time.time() - start_time))
        return polysel1_bestpolys

    def generate_polysel2_progparams(self, params, polysel1_bestpolys):
        I = params["I"]
        alim = params["alim"]
        rlim = params["rlim"]
        N = params["N"]
        outputdir = params["outputdir"]

        progparams = params.copy()
        progparams.setdefault("area", 2. ** (2 * I - 1) * alim)
        progparams.setdefault("Bf", float(alim))
        progparams.setdefault("Bg", float(alim))
        progparams.pop("I", None)
        progparams.pop("alim", None)
        progparams.pop("rlim", None)
        progparams.pop("N", None)
        progparams.pop("batch_size", None)
        progparams.pop("outputdir", None)

        idx = 0
        for poly in polysel1_bestpolys:
            polyfile = str(os.path.join(self.workdir, "%s.polyselect2.raw_%d" % (self.name, idx)))
            outputfile = str(os.path.join(outputdir, "%s.polyselect2.opt_%d" % (self.name, idx)))
            poly.create_file(polyfile)
            program = cadoprograms.PolyselectRopt(inputpolys=polyfile, stdout=outputfile, **progparams)
            idx += 1
            task_command = program.make_command_line()
            cmd_logger.debug(task_command)
            yield task_command, outputfile

    def submit_polysel2_batch(self, generator, params):
        batch_size = params["batch_size"]
        batch_file = str(os.path.join(params["outputdir"], "polysel2.sh"))
        task_commands = []
        for i in range(batch_size):
            try:
                task_commands.append(next(generator))
            except StopIteration:
                break
        logger.info("Submitting %d additional polysel2 tasks", len(task_commands))

        outputfiles = []
        workdir = params.get('workdir')

        i = 0
        num_submitted = 0
        for command, outputfile in task_commands:
            with open(batch_file, 'w') as f:
                outputfiles.append(outputfile)
                f.write("#!/bin/sh\n")
                f.write("#SBATCH -p factor\n")
                f.write("#SBATCH -J %s\n" % outputfile.split('/')[-1])
                f.write("#SBATCH -n 1\n")
                f.write("#SBATCH -c 2\n")
                f.write("#SBATCH -s\n")
                f.write("#SBATCH --requeue\n")
                f.write("#SBATCH --output=/dev/null\n")
                f.write("srun -s -c 2 -n 1 -p factor -J polysel2 " + command + " 2>&1 &\n")
                f.write("wait\n")
            os.chmod(batch_file, 0o755)
            # TODO: use the stdout of sbatch to determine job IDs, and only continue when all of the jobs have finished.
            stdout = utils.run_command("sbatch " + batch_file)

            # rate-limit job submission so slurm is not overwhelmed
            if i >= 100:
                logger.debug("Submitted %d/%d polysel2 jobs", num_submitted, len(task_commands))
                #time.sleep(1)
                i = 0
            num_submitted += 1
            i += 1

        while (True):
            jobs_out = utils.run_command("squeue -t PENDING,RUNNING,COMPLETING").strip().split('\n')
            if len(jobs_out) - 1 <= 0:
                break
            logger.info("Number of queued jobs: %d", len(jobs_out) - 1)
            sleep(10)

        polys = []
        for filen in outputfiles:
            try:
                with open(filen) as f:
                    polys += list(utils.parse_poly(f.read()))
            except Exception as e:
                logger.error("%s. Is workdir NFS-shared?", e)
                raise
        return polys

    def run_polysel2(self, parameters, polysel1_bestpolys):
        # get parameters for polysel2
        polysel2_paths = ["tasks", "polyselect", "polyselect2", "polyselect_ropt"]
        polysel2_program = cadoprograms.PolyselectRopt
        polysel2_keys = polysel2_program.get_accepted_keys()
        polysel2_keys.update({"batch_size": int, "N": int, "I": int, "alim": int, "rlim": int})
        polysel2_params = parameters.myparams(polysel2_keys, polysel2_paths)
        polysel2_params.update({"outputdir": str(os.path.join(self.workdir, self.name + ".upload"))})
        if not os.path.exists(polysel2_params["outputdir"]):
            logger.info("Creating directory for polysel2 output files %s", polysel2_params["outputdir"])
            os.makedirs(polysel2_params["outputdir"])

        logger.info("Starting polysel2")
        start_time = time.time()

        poly_path = os.path.join(self.workdir, "polysel2_input_%d") % (time.time())
        generator = self.generate_polysel2_progparams(polysel2_params, polysel1_bestpolys)

        polys = []
        while True:
            new_polys = self.submit_polysel2_batch(generator, polysel2_params)
            if len(new_polys) == 0:
                break
            polys += new_polys
        
        polysel2_bestpoly = max(polys, key=lambda s: s.MurphyE)
        logger.debug("Polysel2 best poly: %s", polysel2_bestpoly)

        poly_file = os.path.join(self.workdir, self.name + ".polyselect2.poly")
        polysel2_bestpoly.create_file(poly_file)
        logger.info("Polysel2: Saving best polynomial at %s", poly_file)
        logger.info("Polysel2 finished in %s", utils.str_time(time.time() - start_time))
        return poly_file
      
    def msieve_convert_poly(self, polyfile):
        # convert polynomial to msieve format
        logger.info("--- Converting polynomial to msieve format ---")
        binary_path = self.parameters.myparams({'execpath': str}, ["tasks", "polyselect"])['execpath'] + "/misc/convert_poly"

        if not os.path.isfile(polyfile):
            logger.error("Could not find the .polyselect2.poly polynomial file.")
            sys.exit(1)
        if not os.path.isfile(binary_path):
            logger.error("Could not find the convert_poly binary in cado_bindir.")
            sys.exit(1)
        cmd = "{binary} -of msieve < {poly} > {fb}".format(binary=binary_path, poly=polyfile, fb=self.msieve_fb_path)
        cmd_logger.debug(cmd)
        ret = subprocess.call(cmd, shell=True)
        if ret != 0:
            logger.error("Could not convert polynomial. Exiting...")
            sys.exit(1)
        logger.info("--- End of converting polynomial to msieve format ---")

    def run(self):        
        # Check if we should import a polynomial
        import_poly = self.parameters.myparams({'import': [None]}, ['tasks', 'polyselect']).get('import')
        polyfile = None
        if import_poly:
            logger.info("Importing polynomial from '%s'", import_poly)
            polyfile =  import_poly
        else:
            polysel1_bestpolys = self.run_polysel1(self.parameters)
            polyfile = self.run_polysel2(self.parameters, polysel1_bestpolys)

        # if we are using msieve for linear algebra, convert the polynomial to the correct format
        self.msieve_convert_poly(polyfile)

        return polyfile
