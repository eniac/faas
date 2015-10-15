import sys
import errno
from collections import defaultdict 
import os
import math
import re
import time
import utils
import cadotask
from cadotask import Duplicates1Task, Duplicates2Task, Polynomials
import cadoprograms
import cadoparams

import logging
logger = logging.getLogger('Filter')
cmd_logger = logging.getLogger('cmd')

def init(params, freerel):
    global parameters
    parameters = params
    global freerel_out
    freerel_out = freerel
    p = parameters.myparams({'workdir': str, 'name': str}, ['tasks'])
    name = p.get('name')
    workdir = p.get('workdir')
    # convenient global variables
    global dup1dir
    dup1dir = os.path.join(workdir, name + '.duplicates1')
    global log_path
    log_path = os.path.join(workdir, name + "_msieve_log_dir")
    global dat_path
    dat_path = os.path.join(workdir, name + ".dat")
    global cyc_path
    cyc_path = dat_path + ".cyc"
    global fb_path
    fb_path = os.path.join(workdir, name + ".fb")
    global purgedfile
    purgedfile = os.path.join(workdir, "purge.purged.gz")

# takes the output from las, and outputs the number of relations found
def parse_slice_counts(stderr):
    """ Takes lines of text and looks for slice counts as printed by dup1 """
    counts = [None] * generate_dup1_command.nr_slices
    for line in stderr.splitlines():
        match = re.match(r'# slice (\d+) received (\d+) relations', line)
        if match:
            (slicenr, nrrels) = map(int, match.groups())
            if not counts[slicenr] is None:
                raise Exception("Received two values for relation count "
                    "in slice %d" % slicenr)
            counts[slicenr] = nrrels
    for (slicenr, nrrels) in enumerate(counts):
        if nrrels is None:
            return 0
            raise Exception("Received no value for relation count in "
                "slice %d" % slicenr)
    return counts

# Parses the output from cadonfs purge, and returns the number of relations still needed
def parse_purge_stderr(stderr, input_nrels):
    # If stderr ends with
    # b'excess < 0.10 * #primes. See -required_excess argument.'
    # then we need more relations from filtering and return False

    not_enough1 = re.compile(r"\(excess / nprimes\) = \d+.?\d* < \d+.?\d*. "
        r"See -required_excess argument.")
    not_enough2 = re.compile(r"number of rows <= number of columns")

    nrels_nprimes = re.compile(r"\s*nrows=(\d+)\sncols=(\d+)\sexcess=(-?\d+)")

    # if you don't find either of the not enough lines, we have required rels
    if not (re.findall(not_enough1, stderr) or re.findall(not_enough2, stderr)):
        return 0
    matches = re.findall(nrels_nprimes, stderr)
    assert len(matches) > 2

    (input_nrels, input_nprimes, input_excess) = map(int, matches[0])
    (nrels, nprimes, excess) = map(int, matches[-1])
    logger.debug("input_excess: %d", input_excess)
    logger.debug("nrels: %d\tnprimes: %d\t excess: %d", nrels, nprimes, excess)
    assert nrels - nprimes == excess

    if excess > 0:
        return 0

    add = math.ceil(input_nrels / nrels * -excess) * 1.1
    add = min(add, input_nrels * 2)
    return int(add)


def generate_dup1_command(newfiles):
    paths = ["tasks", "filter", "duplicates1", "dup1"]
    program_to_run = cadoprograms.Duplicates1

    progparams = parameters.myparams(program_to_run.get_accepted_keys(), paths)
    progparams.pop("prefix", None)
    progparams.pop("out", None)

    prefix = "dup1.%d" % (time.time())

    if generate_dup1_command.nr_slices is None:  # first run
        generate_dup1_command.nr_slices = 2  # default to 2**1
        if "nslices_log" in progparams:
            generate_dup1_command.nr_slices = 2 ** progparams["nslices_log"]
        for i in range(0, generate_dup1_command.nr_slices):
            try:
                os.makedirs(os.path.join(dup1dir, str(i)))
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise

    # might want to improve this and maintain status of all files 
    # (in case a dup1 run fails) so that a later run can try again
    if len(newfiles) <= 10:
        program = cadoprograms.Duplicates1(*newfiles, prefix=prefix, out=dup1dir, **progparams)
    else:
        filelist = utils.write_list_to_file(newfiles, dup1dir)
        program = cadoprograms.Duplicates1(filelist=filelist, prefix=prefix, out=dup1dir, **progparams)

    command = program.make_command_line()
    cmd_logger.debug(command)
    return command
# hacky way to persist this without a class
generate_dup1_command.nr_slices = None


# Returns a list of tuples in the format (slice_filename: number_of_relations_in_slice)
def run_dup1(newfiles):
    cmd = generate_dup1_command(newfiles)
    stdout, stderr = utils.run_command(cmd, include_stdout=True, include_stderr=True, logger=logger)
    files = Duplicates1Task.parse_output_files(stderr)  # filenames are in a dict{name:slice_idx}
    sorted_files = sorted(files.keys(), key=files.get)
    counts = parse_slice_counts(stderr)
    if not counts:
        return -1

    return list(zip(sorted_files, counts))


def generate_dup2_commands(newfiles):
    renumberfilename = freerel_out["renumberfilename"]
    paths = ["tasks", "filter", "duplicates2", "dup2"]
    program_to_run = cadoprograms.Duplicates2
    progparams = parameters.myparams(program_to_run.get_accepted_keys(), paths)
    progparams.pop("renumber", None)

    commands = []

    for i, (filename, rels) in enumerate(newfiles):
        if filename in generate_dup2_commands.slice_files[i]:
            continue  # nothing to do

        logger.info("Dup2: Processing slice %d", i)
        generate_dup2_commands.slice_rels[i] += rels
        rels = generate_dup2_commands.slice_rels[i] 
        generate_dup2_commands.slice_files[i].add(filename)
        files = list(generate_dup2_commands.slice_files[i])

        if len(files) <= 10:
            program = cadoprograms.Duplicates2(*files, rel_count=rels, renumber=renumberfilename, **progparams)
        else:
            filelist = utils.write_list_to_file(files, dup1dir)
            program = cadoprograms.Duplicates2(filelist=filelist, rel_count=rels, renumber=renumberfilename, **progparams)
        command = program.make_command_line()
        cmd_logger.debug(command)
        commands.append(command)
    return commands

# hacky way to persist this without a class
generate_dup2_commands.slice_files = defaultdict(set)
generate_dup2_commands.slice_rels = defaultdict(int)


def run_dup2(newfiles):
    task_commands = generate_dup2_commands(newfiles)

    rel_counts = []
    # could parallelize this part too, but # of slices is usually low (2)
    for cmd in task_commands:
        stdout, stderr = utils.run_command(cmd, include_stdout=True, include_stderr=True, logger=logger)
        try:
            rel_counts.append(Duplicates2Task.parse_remaining(stderr.splitlines()))
        except:
            return -1
    return rel_counts

def generate_purge_command(dup2_out):
    paths = ["tasks", "filter", "purge", "purge"]
    program_to_run = cadoprograms.Purge

    progparams = parameters.myparams(program_to_run.get_accepted_keys(), paths)

    nfree = freerel_out["nfree"]
    nunique = sum(dup2_out)
    input_nrels = nfree + nunique
    nprimes = freerel_out["nprimes"]

    minindex = int(progparams.get("col_minindex", -1))
    if minindex == -1:
        minindex = int(nprimes / 20.0)
        # For small cases, we want to avoid degenerated cases, so let's
        # keep most of the ideals: memory is not an issue in that case.
        if (minindex < 10000):
            minindex = 500
        progparams.setdefault("col_minindex", minindex)
    keep = progparams.pop("keep", None)
    
    relsdelfile = None  # not supporting dlp yet

    files = [freerel_out["freerelfilename"]]
    for i in range(generate_dup1_command.nr_slices):
        files += list(generate_dup2_commands.slice_files[i])

    if len(files) <= 10:
        program = cadoprograms.Purge(*files, nrels=input_nrels, out=purgedfile, outdel=relsdelfile, keep=keep, nprimes=nprimes, **progparams)
    else:
        filelist = utils.write_list_to_file(files, dup1dir)
        program = cadoprograms.Purge(nrels=input_nrels, out=purgedfile, outdel=relsdelfile, keep=keep, nprimes=nprimes, filelist=filelist, **progparams)

    command = program.make_command_line()
    cmd_logger.debug(command)
    return command, purgedfile

# Calculates how many relations are needed to generate the final matrix
def run_purge(dup2_out):

    cmd, purgedfile = generate_purge_command(dup2_out)
    stdout, stderr = utils.run_command(cmd, include_stdout=True, include_stderr=True, logger=logger)
    nfree = freerel_out["nfree"]
    nunique = sum(dup2_out)
    input_nrels = nfree + nunique

    # parse_purge_stderr(stdout) is not a bug, purge currently outputs everything to stdout
    additional = parse_purge_stderr(stdout, input_nrels)

    return additional

def run_msieve_filter():
    logger.info("--- Beginning MSieve de-duplication and filtering ---")
    params = parameters.myparams({
        "N": int,
        "execpath": str,
        "msievedir": str,
        "target_density": 70,
        "threads": int,
    }, ['tasks', 'msieve', 'filter'])

    number = str(params['N'])
    execpath = params['execpath']
    msievedir = params['msievedir']
    target_density = params['target_density']
    threads = params['threads']

    filter_cmd = msievedir + "/msieve"
    filter_cmd += " -nf " + fb_path
    filter_cmd += " -nc1 \"target_density={td}\"".format(td=target_density)
    filter_cmd += " -v -t " + str(threads)
    filter_cmd += " -l " + log_path
    filter_cmd += " -s " + dat_path
    filter_cmd += " " + number

    cmd_logger.debug(filter_cmd)
    stdout, stderr, ret = utils.run_command(filter_cmd, include_stdout=True, include_stderr=True, include_returncode=True, logger=logger)

    to_return = 0
    if ret != 0:
        logger.error("Received error code " + str(ret) + " from Msieve filtering. Exiting...")
        sys.exit(1)
    elif not os.path.isfile(cyc_path):

        second_to_last_line = stdout.strip().split('\n')[-2]
        match = re.search('filtering wants (\d+) more relations', second_to_last_line)
        if match:
            logger.warning(".cyc file not found. Likely there are too few relations. Sieve more!")
            rels_additional = int(match.group(1))
            to_return = rels_additional
        elif re.search('factor', second_to_last_line):
            # check if msieve just completed the factorization using trial division alone
            logger.info("Filtering completed factorization using trial division")
            to_return = -2
        elif re.search('found (\d+) cycles, need (\d+)', second_to_last_line):
            logger.info('Not enough cycles for linear algebra. Requesting 10% more relations.')
            to_return = -1

    logger.info("--- End of MSieve de-duplication and filtering ---")
    return to_return

def run(relation_files):
    if not relation_files:
        return -1

    if run.do_cado_filtering:
    
        logger.info("--- Beginning Cado filtering ---")
        logger.info('Dup1: Partitioning new files across slices')
        dup1_out = run_dup1(relation_files)
        if dup1_out == -1:
            return -1
        for fn, count in dup1_out:
            if count == 0:
                return -1
        logger.debug("Dup1 out: %s", dup1_out)

        logger.info('Dup2: Removing duplicates')
        dup2_out = run_dup2(dup1_out)
        if dup2_out == -1:
            return -1, None
        logger.debug("Dup2 out: %s", dup2_out)
        # NOTE: purge (and freerel) do not support GaloisFilter
        new_rels_wanted = run_purge(dup2_out)

        # request a minimum of 10k new relations if we are requesting more
        if new_rels_wanted > 0:
            new_rels_wanted = max(new_rels_wanted, 10000)

        logger.info('Purge: %d additional rels requested', new_rels_wanted)
        if new_rels_wanted == 0:
            run.do_cado_filtering = False
            run.do_msieve_filtering = True
            logger.info("--- Finished Cado filtering ---")

    if run.do_msieve_filtering:
        new_rels_wanted = run_msieve_filter()

    return new_rels_wanted

# NOTE: Enable cado filtering by uncommenting the following two lines, and commenting out the two lines after that.
#run.do_cado_filtering = True
#run.do_msieve_filtering = False
run.do_cado_filtering = False
run.do_msieve_filtering = True
