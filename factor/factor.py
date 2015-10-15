import argparse
import os
import time
import utils
import polysel
import sieving
import linalg
import sqrt
import math
import cadoparams
import cadoprograms

import re

# signal handling
import signal
import sys

import logging
import types

logger = logging.getLogger('Factor')
cmd_logger = logging.getLogger('cmd')

# if the progarm is terminated by Ctrl +C, send a kill signal to running slurm jobs
def signal_handler(signal, frame):
    logger.warning("Interrupt received. Terminating...")
    utils.run_command("scancel -p factor")
    sys.exit(1)

def setup_logging(workdir, name):
    '''
    LOGGING LEVELS:
        CRITICAL    50
        ERROR   40
        WARNING 30
        INFO    20
        DEBUG   10
        NOTSET  0
    '''

    # Formats for logging
    FILE_FORMAT='PID%(process)d %(asctime)-15s %(name)s:%(levelname)s: %(message)s'
    CONSOLE_FORMAT='%(name)s:%(levelname)s: %(message)s'
    CMD_FORMAT='#PID%(process)d %(asctime)-15s\n%(message)s\n'
    
    # Set the log level for messages in the log file
    file_handler = logging.FileHandler(os.path.join(workdir, name + '.log'))
    file_handler.setLevel(logging.NOTSET)
    file_formatter = logging.Formatter(FILE_FORMAT)
    file_handler.setFormatter(file_formatter)
    logging.getLogger('').addHandler(file_handler)

    # Set the log level for messages printed to the console
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.NOTSET)
    console_formatter = logging.Formatter(CONSOLE_FORMAT)
    console_handler.setFormatter(console_formatter)
    logging.getLogger('').addHandler(console_handler)

    # Log all commands to file with the cmd logger
    cmd_logger = logging.getLogger('cmd')
    cmd_handler = logging.FileHandler(os.path.join(workdir, name + '.cmd'))
    cmd_handler.setLevel(logging.NOTSET)
    cmd_formatter = logging.Formatter(CMD_FORMAT)
    cmd_handler.setFormatter(cmd_formatter)
    cmd_logger.addHandler(cmd_handler)
    cmd_logger.propagate = False

    # Also log cmds to the log file and to stdout
    #cmd_logger.addHandler(file_handler)
    #cmd_logger.addHandler(console_handler)

# Estimate the number of relations that will be required to pass filtering
# based on previous experiments.
def estimate_rels_wanted(parameters):
    p = parameters.myparams({'lpba':int,'lpbr':int,'target_density':int, 'N':int}, ['tasks', 'msieve'])
    d = {}
    td = p['target_density']
    n = len(str(p['N']))
    if n == 155 or n == 154: # 512-bit numbers can have either 154 or 155 digits
        if p['lpba'] == 28 and p['lpbr'] == 28:
            d = {
                    70: 28200000,
                    80: 28200000,
                    90: 28300000,
                    100: 28400000,
                    110: 28900000,
                    120: 30300000,
                    130: 33100000,
                 }

        elif p['lpba'] == 29 and p['lpbr'] == 29:
            d = {
                    70: 45300000,
                    80: 45400000,
                    90: 45500000,
                    100: 45900000,
                    110: 46800000,
                    120: 48600000,
                    130: 51500000,
                    140: 57800000,
                }
    if n == 120:
        if p['lpba'] == 26 and p['lpbr'] == 26:
            d = {
                    70: 7450000, 
                    80: 7450000,
                    90: 7450000,
                    100: 7450000,
                    110: 9500000,
                    120: 13500000,

                }
    
    # round up to closest target density to get the initial number of relations
    num_rels = 0
    for i in sorted(d.keys()):
        num_rels = d[i]
        if i > td:
            break
    return num_rels

# Add new parameters that may not have been present in the parameters file.
# NOTE: This currently will overwrite parameters from the parameter file.
def set_static_parameters(parameters):
    name = parameters.myparams({"name": str}, ['tasks'])["name"]
    workdir = parameters.myparams({"workdir": str}, ['tasks'])["workdir"]

    if not parameters.myparams({'target_density': [int]}, ['tasks', 'msieve']):
        parameters.readparams(['tasks.msieve.target_density=70'])

    # If rels_wanted is not set, and we have an idea of the minimum number of relations required for this target density, set rels_wanted
    if parameters.myparams({'rels_wanted': 0}, ['tasks', 'sieve'])['rels_wanted'] == 0:
        rels_wanted = estimate_rels_wanted(parameters)
        parameters.readparams(['tasks.sieve.rels_wanted=%d' % rels_wanted])
        logger.info('rels_wanted unspecified. Setting rels_wanted based on target_density, lpba, and lpbr (%d)', rels_wanted)

    cores_info = {k: v for k, v in (x.split('=') for x in utils.run_command("slurmd -C").split())}
    spm = int(cores_info['Boards']) * int(cores_info['SocketsPerBoard'])
    cps = int(cores_info['CoresPerSocket'])
    tpc = int(cores_info['ThreadsPerCore'])

    mpi_rows, mpi_cols = [int(x) for x in parameters.myparams({"mpi": str}, ['tasks', 'msieve', 'linalg'])['mpi'].split("x")]
    mpi_hosts = parameters.myparams({'hosts': str}, ['tasks', 'msieve', 'linalg'])['hosts'].split(',')

    # checking the num of mpi jobs does not exceed total core count
    if (mpi_cols * mpi_rows > cps * spm * len(mpi_hosts)):
        logger.error("Cannot have more jobs per machine than physical cores. Exiting...")
        sys.exit(1)

    new_params = [
            "tasks.msieve.rankfile = %s" % os.path.join(workdir, name + ".rankfile"),
            "tasks.msieve.fb_path = %s" % os.path.join(workdir, name + ".fb"),
            "tasks.msieve.dat_path = %s" % os.path.join(workdir, name + ".dat"),
            "tasks.msieve.log_path = %s" % os.path.join(workdir, name + ".msieve.log"),
            "tasks.msieve.linalg.phys_socks_per_machine = %d" % spm,
            "tasks.msieve.linalg.phys_core_per_sock = %d" % cps,
            "tasks.msieve.linalg.threads_per_core = %d" % tpc,
            "tasks.filter.threads = %d" % max(int(tpc*cps*spm),36),
            "tasks.msieve.filter.threads = %d" % max(int(tpc*cps*spm),36),
            "tasks.msieve.linalg.threads = %d" % int((math.floor(tpc*cps*spm / math.ceil(mpi_rows*mpi_cols) / len(mpi_hosts)))),
            "tasks.msieve.sqrt.threads = %d" % max(int(tpc*cps*spm),36),
            ]

    parameters.readparams(new_params)
    return parameters

# Set parameters based on the cluster setup (e.g., total number of available cores).
def set_dynamic_parameters(parameters):
    # adjust parameters to fit the cluster that factorization is being run on.
    # if unspecified, set the slurm batch size
    if parameters.myparams({'batch_size': 0}, ['tasks'])['batch_size'] == 0:
        info = utils.run_command("sinfo -p factor -o %C --states=ALLOC,ALLOCATED,COMP,COMPLETING,IDLE").split('\n')[1].split('/')
        alloc_cores = int(info[0])
        free_cores = int(info[1])
        num_cores = alloc_cores + free_cores
        parameters.readparams(['tasks.batch_size=%d' % num_cores])
        logger.info("Batch size unspecified. Setting batch_size = #cores (%d)", num_cores)
        # Limit the polyselect batch size, since more cores don't seem to matter at this point
        parameters.readparams(['tasks.polyselect.batch_size=%d' % min(1000,num_cores)])

    # Set the adrange and nrkeep parameters to match batch_size
    admax = parameters.myparams({'admax': int}, ['tasks', 'polyselect'])['admax']
    polysel_batch = parameters.myparams({'batch_size': int}, ['tasks', 'polyselect'])['batch_size']
    adrange = int(math.ceil(int(admax)/polysel_batch))
    parameters.readparams(['tasks.polyselect.adrange=%d' % adrange, 'tasks.polyselect.nrkeep=%d' % polysel_batch])
    logger.info("Setting adrange = admax/#cores (%d), and nrkeep = #cores (%d)", adrange, num_cores)

    return parameters

# Check that all of the required parameters are present, and return a dictionary of { stage: stage_params }
def check_parameters(parameters):
    params = {}
    # check general parameters
    stage = 'start'
    stage_paths = ['tasks']
    stage_keys = {"name": str, "workdir": str, "N": int, "execpath": str}
    stage_params = parameters.myparams(stage_keys, stage_paths)
    params[stage] = stage_params

    # check polysel parameters
    stage = 'polysel'
    stage_paths = ['tasks', 'polyselect']
    stage_keys = {'import': [str]}
    stage_params = parameters.myparams(stage_keys, stage_paths)
    params[stage] = stage_params

    # check polysel1 parameters
    stage = 'polysel1'
    stage_paths = ["tasks", "polyselect", "polyselect1", "polyselect2l"]
    stage_program = cadoprograms.Polyselect2l
    stage_keys = stage_program.get_accepted_keys()
    stage_params = parameters.myparams(stage_keys, stage_paths)
    stage_params.update({"outputdir": str})
    params[stage] = stage_params

    # check polysel2 parameters
    stage = 'polysel2'
    stage_paths = ["tasks", "polyselect", "polyselect2", "polyselect_ropt"]
    stage_program = cadoprograms.PolyselectRopt
    stage_keys = stage_program.get_accepted_keys()
    stage_keys.update({"I": int, "alim": int, "rlim": int})
    stage_params = parameters.myparams(stage_keys, stage_paths)
    stage_params.update({"outputdir": str})
    params[stage] = stage_params

    # check sieve parameters
    stage = 'sieve'
    stage_paths = ["tasks", "sieve"]
    stage_keys = {"import": [None], "alim": int, "rlim": int, "qrange": int, "lpba": int, "qmin": [int], "execpath": str, "rels_wanted": int}
    stage_params = parameters.myparams(stage_keys, stage_paths)
    params[stage] = stage_params

    # check las parameters
    stage = 'las'
    stage_paths = ["tasks", "sieve", "las"]
    stage_program = cadoprograms.Las
    stage_keys = stage_program.get_accepted_keys()
    stage_params = parameters.myparams(stage_keys, stage_paths)
    params[stage] = stage_params

    # check factorbase parameters
    stage = 'factorbase'
    stage_paths = ["tasks", "sieve", "factorbase", "makefb"]
    stage_program = cadoprograms.MakeFB
    stage_keys = stage_program.get_accepted_keys()
    stage_keys.update({"gzip": True, "I": int, "rlim": int, "alim": int})
    stage_params = parameters.myparams(stage_keys, stage_paths)
    params[stage] = stage_params

    # check freerel parameters 
    stage = 'freerel'
    stage_paths = ["tasks", "sieve", "freerel"]
    stage_program = cadoprograms.FreeRel
    stage_keys = stage_program.get_accepted_keys()
    stage_params = parameters.myparams(stage_keys, stage_paths)
    params[stage] = stage_params

    # check dup1 parameters
    stage = 'dup1'
    stage_paths = ["tasks", "filter", "duplicates1", "dup1"]
    stage_program = cadoprograms.Duplicates1
    stage_keys = stage_program.get_accepted_keys()
    stage_params = parameters.myparams(stage_keys, stage_paths)
    params[stage] = stage_params

    # check dup2 parameters
    stage = 'dup2'
    stage_paths = ["tasks", "filter", "duplicates2", "dup2"]
    stage_program = cadoprograms.Duplicates1
    stage_keys = stage_program.get_accepted_keys()
    stage_params = parameters.myparams(stage_keys, stage_paths)
    params[stage] = stage_params

    # check purge parameters
    stage = 'purge'
    stage_paths = ["tasks", "filter", "purge", "purge"]
    stage_program = cadoprograms.Purge
    stage_keys = stage_program.get_accepted_keys()
    stage_params = parameters.myparams(stage_keys, stage_paths)
    params[stage] = stage_params

    # check msieve parameters
    stage = 'msieve'
    stage_paths = ['tasks', 'msieve']
    stage_keys = {'target_density': int, 'msievedir': str}
    stage_params = parameters.myparams(stage_keys, stage_paths)
    params[stage] = stage_params

    # check msieve filter parameters
    stage = 'msieve_filter'
    stage_paths = ['tasks', 'msieve', 'filter']
    stage_keys = {'threads': int, 'dat_path': str, 'log_path': str} 
    stage_params = parameters.myparams(stage_keys, stage_paths)
    params[stage] = stage_params

    # check msieve linalg parameters
    stage = 'msieve_linalg'
    stage_paths = ['tasks', 'msieve', 'linalg']
    stage_keys = {"mpi": str, "hosts": str, "phys_socks_per_machine": int, "phys_core_per_sock": int, "threads_per_core": int}
    stage_params = parameters.myparams(stage_keys, stage_paths)
    params[stage] = stage_params

    # check msieve sqrt parameters
    stage = 'msieve_sqrt'
    stage_paths = ['tasks', 'msieve', 'sqrt']
    stage_keys = {"threads": int, "log_path": str, "fb_path": str, "dat_path": str}
    stage_params = parameters.myparams(stage_keys, stage_paths)
    params[stage] = stage_params

    return params

def checkpoint_stage(stage, result):
    checkpoint = utils.get_checkpoint().get(stage, {})
    checkpoint.update({'result': result})
    utils.update_checkpoint({
        'stage': stages.index(stage),
        stage: checkpoint,
        })

def load_stage(stage):
    return utils.get_checkpoint().get(stage, {}).get('result')

def stage_required(stage):
    # check that 'stage-1' has been checkpointed
    stage_index = stages.index(stage)
    required = (utils.get_checkpoint()['stage'] < stage_index)

    if stage_required.manual_stages:
        if stage_index > max(stage_required.manual_stages):
            sys.exit(0)
        elif stage_index in stage_required.manual_stages:
            return True
        else:
            return False
    return required
stage_required.manual_stages = []

def set_manual_stages(params):
    start_stage_idx = len(stages) - 1
    checkpoint = utils.get_checkpoint()['params']

    # For each stage, check if any of the parameters have changed since the last run
    changed_stages = set()
    for stage,stage_params in params.items():
        diff = {}
        for key,value in stage_params.items():
            if checkpoint.get(stage) and checkpoint.get(stage).get(key) != value:
                # the value differs from the checkpoint
                diff[key] = (value, checkpoint.get(stage).get(key))

        if diff:
            logger.info("parameters differ from checkpoint in stage %s: %s", stage, diff)
            changed_stages.add(stage)

    # For stage, mark is as changed if any of its dependencies changed
    for stage in stages:
        for dependency in dependencies.get(stage, []):
            if stages[dependency] in changed_stages:
                changed_stages.add(stage)
    
    # Finally, add the required stages to the manual stages
    if changed_stages:
        start_stage_idx = min([ stages.index(stage) for stage in changed_stages ])
    stage_required.manual_stages = range(start_stage_idx, len(stages))

# Each of the checkpointed stages
stages = ['start','polysel', 'polysel1', 'polysel2', 'sieve', 'las', 'factorbase', 'freerel', 'dup1', 'dup2', 'purge', 'cado_filter', 'msieve', 'msieve_filter', 'linalg', 'msieve_linalg', 'sqrt', 'msieve_sqrt', 'complete']

# If one of the dependencies for a particular stage is run again, then that stage should also be run again.
dependencies = {
        'start': range(stages.index('start'), stages.index('start')),
        'polysel': range(stages.index('start'), stages.index('polysel2')), 
        'sieve': range(stages.index('start'), stages.index('msieve_filter')),
        'linalg': range(stages.index('start'), stages.index('msieve_linalg')), 
        'sqrt': range(stages.index('start'), stages.index('msieve_sqrt')),
        'complete': range(stages.index('start'), stages.index('complete')),
        }
            
def do_polysel(parameters):
    name = parameters.myparams({"name": str}, ['tasks'])["name"]
    workdir = parameters.myparams({"workdir": str}, ['tasks'])["workdir"]
    result = {}
    if stage_required('polysel'):
        poly_start = time.time()
        polysel_offset = utils.slurm_cpu_time_start()
        p = polysel.Polysel(parameters)
        result['polyfile'] = p.run()
        poly_finish = time.time()
        result['duration'] = poly_finish - poly_start
        logger.info("\tPolysel in %s", utils.str_time(result['duration']))

        if parameters.myparams({"collect_cputime": [str]}, ['tasks', 'polyselect']).get("collect_cputime"):
            logger.info("--- Collecting polysel cumulative CPUTime ---")
            polysel_slurm_outfile = str(os.path.join(workdir, name + ".slurm_polysel"))
            result['cputime'] = utils.slurm_cpu_time_end(polysel_offset, polysel_slurm_outfile)
            logger.info("\tPolysel cumulative CPUTime %s" % utils.str_time(result['cputime']))
        else:
            result['cputime'] = 0
        
        checkpoint_stage('polysel', result)
    else:
        result = load_stage('polysel')

    return result

def do_sieve(parameters, polysel_result):
    result = {}
    name = parameters.myparams({"name": str}, ['tasks'])["name"]
    workdir = parameters.myparams({"workdir": str}, ['tasks'])["workdir"]
    # msieve_filter is the last checkpointed stage in this block
    if stage_required('sieve'):
        sieve_start = time.time()
        sieving_offset = utils.slurm_cpu_time_start()
        s = sieving.Sieve(parameters, polysel_result['polyfile'])
        relation_files = s.run()
        result['relation_files'] = relation_files
        sieve_finish = time.time()
        result['duration'] = sieve_finish - sieve_start
        logger.info("\tSieving in %s", utils.str_time(result['duration']))

        if parameters.myparams({"collect_cputime": [str]}, ['tasks', 'sieve']).get("collect_cputime"):
            logger.info("--- Collecting sieving cumulative CPUTime ---")
            sieving_slurm_outfile = str(os.path.join(workdir, name + ".slurm_sieving"))
            result['cputime'] = utils.slurm_cpu_time_end(sieving_offset, sieving_slurm_outfile)
            logger.info("\tSieving cumulative CPUTime %s" % utils.str_time(result['cputime']))
        else:
            result['cputime'] = 0

        relation_files_file = str(os.path.join(workdir, name + '.relation_files'))
        logger.info("Saving used relation files from this run to %s", relation_files_file)
        with open(relation_files_file, 'w') as f:
            for line in relation_files:
                f.write(line+'\n')

        if s.completed_factorization:
            logger.info("Completed factorization using trial division. Check log file for factors.")
            utils.update_checkpoint({
                'trail_division': True})

        checkpoint_stage('sieve', result)
    else:
        result = load_stage('sieve')

    post_sieve = parameters.myparams({'post_sieve': None}, ['commands']).get('post_sieve')
    if post_sieve != None:
        logger.info('Post-sieve command %s', post_sieve)
        utils.run_command(post_sieve, logger=logger)
    
    return result

def do_linalg(parameters, sieve_result):
    result = {}
    # Check the rare case that factorization was completed in msieve's filtering
    if utils.get_checkpoint().get('trial_division') != None:
        result['duration'] = 0
    elif stage_required('linalg'):
        linalg_start = time.time()
        linalg.run(parameters)
        linalg_finish = time.time()
        result['duration'] = linalg_finish - linalg_start
        logger.info("\tLinalg in %s", utils.str_time(result['duration']))

        checkpoint_stage('linalg', result)
    else:
        result = load_stage('linalg')

    post_linalg = parameters.myparams({'post_linalg': None}, ['commands']).get('post_linalg')
    if post_linalg != None:
        logger.info('Post-linalg command %s', post_linalg)
        utils.run_command(post_linalg, logger=logger)

    return result

def do_sqrt(parameters, linalg_result):
    result = {}
    if stage_required('sqrt'):

        sqrt_start = time.time()
        result['factors'] = sqrt.run(parameters)
        sqrt_finish = time.time()
        result['duration'] = sqrt_finish - sqrt_start
        logger.info("\tSqrt in %s", utils.str_time(result['duration']))

        checkpoint_stage('sqrt', result)
    else:
        result = load_stage('sqrt')

    post_sqrt = parameters.myparams({'post_sqrt': None}, ['commands']).get('post_sqrt')
    if post_sqrt != None:
        logger.info('Post-sqrt command %s', post_sqrt)
        utils.run_command(post_sqrt, logger=logger)

    return result

def main():
    
    signal.signal(signal.SIGINT, signal_handler)

    parser = argparse.ArgumentParser(description="Integer Factorization with "
                                         "the Number Field Sieve")
    parser.add_argument("parameters", help="A file with the parameters to use")
    parser.add_argument("options", metavar="OPTION", help="An option as in "
                            "parameter file (format: key=value)", nargs="*")
    parser.add_argument('--resume','-r', help="checkpoint file to resume from")
    parser.add_argument('--stage','-s', action='append', help="stage to complete ('start','polysel','sieving','linalg','complete'), add + to run all subsequent stages")
    
    
    args = parser.parse_args()
    parameters = utils.get_params(args.parameters, args.options)

    name = parameters.myparams({"name": str}, ['tasks'])["name"]
    workdir = parameters.myparams({"workdir": str}, ['tasks'])["workdir"]

    if not os.path.exists(workdir):
        logger.info("Creating work directory %s", workdir)
        os.makedirs(workdir)

    setup_logging(workdir, name)

    # Load or create initial checkpoint
    checkpoint_file = args.resume
    if not checkpoint_file:
        checkpoint_file = os.path.join(workdir, "checkpoint.dat")
    utils.init_checkpoint(checkpoint_file)

    # set parameters that are unlikely to change from run to run, such as filenames and directories
    parameters = set_static_parameters(parameters)

    # check that all required parameters are present
    params = check_parameters(parameters)
    utils.update_checkpoint({'params': params})

    # set parameters that will likely change from run to run
    parameters = set_dynamic_parameters(parameters)

    # Write a snapshot of the parameters to a file
    snapshot_filename = "%s/%s.parameters_snapshot" % (workdir, name)
    with open(snapshot_filename, "w") as snapshot_file:
        logger.debug("Writing parameter snapshot to %s", snapshot_filename)
        snapshot_file.write(str(parameters))
        snapshot_file.write("\n")

    start_time = time.time()

    # For each checkpointed stage, check if the stage should be run again.
    # A stage should be run again under the following circumstances:
    #   - The user manually requested to run the stage
    #   - No checkpoint exists for the stage
    #   - A stage on which this stage depends will be re-run
    #   - Parameters on which the stage depends have been changed since the last run
    if args.stage:
        for stage in args.stage:
            if stage.endswith('+'):
                stage = stage[:-1]
                if stage not in stages:
                    continue
                stage_required.manual_stages = range(stages.index(stage), len(stages))
                break
            if stage not in stages:
                args.stage.pop(stage)
            else:
                stage_required.manual_stages.append(stages.index(stage))
    else:
        # since no stage were specified to run manually, choose the first stage based on the checkpoint file
        stage_required.manual_stages = set_manual_stages(params)

    # Run polynomial selection
    polysel_result = do_polysel(parameters)

    # Run sieving
    sieve_result = do_sieve(parameters, polysel_result)

    # Run linalg
    linalg_result = do_linalg(parameters, sieve_result)

    # Run square root
    sqrt_result = do_sqrt(parameters, linalg_result)

    factoring_duration = polysel_result['duration'] + sieve_result['duration'] + linalg_result['duration'] + sqrt_result['duration']
    logger.info('Factoring completed in %s', utils.str_time(factoring_duration))
    logger.info('\tPolysel in real/cpu %s/%s', utils.str_time(polysel_result['duration']), utils.str_time(polysel_result['cputime']))
    logger.info("\tSieving in real/cpu %s/%s", utils.str_time(sieve_result['duration']), utils.str_time(sieve_result['cputime']))
    logger.info("\tLinalg in %s", utils.str_time(linalg_result['duration']))
    logger.info("\tSqrt in %s", utils.str_time(sqrt_result['duration']))
    logger.info("\tFactors %s", ','.join(sqrt_result['factors']))

    post_factor = parameters.myparams({'post_factor': None}, ['commands']).get('post_factor')
    if post_factor != None:
        logger.info('Post-factor command %s', post_factor)
        utils.run_command(post_factor, logger=logger)

if __name__ == "__main__":
    main()
