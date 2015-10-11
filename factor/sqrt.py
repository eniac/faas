import re
import cadoparams
import utils
import os
import sys

import logging
logger = logging.getLogger('Sqrt')
cmd_logger = logging.getLogger('cmd')

def generate_sqrt_command(sqrt_params):
    sqrt_cmd = os.path.join(sqrt_params['msievedir'], 'msieve')
    sqrt_cmd += " -nf " + sqrt_params['fb_path']
    sqrt_cmd += " -nc3 \"target_density={td}\"".format(td=sqrt_params['target_density'])
    sqrt_cmd += " -v "
    sqrt_cmd += " -t " + str(sqrt_params['threads'])
    sqrt_cmd += " -l " + sqrt_params['log_path']
    sqrt_cmd += " -s " + sqrt_params['dat_path']
    sqrt_cmd += " " + str(sqrt_params['N'])
    return sqrt_cmd

def run(parameters):
    # check msieve sqrt parameters
    sqrt_paths = ['tasks', 'msieve', 'sqrt']
    sqrt_keys = {
            "msievedir": str,
            "threads": int, 
            "target_density": int, 
            "fb_path": str, 
            "dat_path": str, 
            "log_path": str, 
            "N": int
            }
    sqrt_params = parameters.myparams(sqrt_keys, sqrt_paths)

    logger.info("--- Beginning MSieve square root phase ---")

    sqrt_cmd = generate_sqrt_command(sqrt_params)

    cmd_logger.debug(sqrt_cmd)
    stdout, stderr, ret = utils.run_command(sqrt_cmd, include_stdout=True, include_stderr=True, include_returncode=True, logger=logger)
    logger.info(stdout)
    logger.info(stdout)
    if ret != 0:
        logger.error("Received error code " + str(ret) + " from Msieve square root. Exiting...")
        sys.exit(1)

    logger.info("--- End of MSieve square root phase ---")
    
    factor_regex = re.compile("^prp\d+ factor: (\d+)$")
    factors = []
    for line in stdout.split('\n'):
        match = factor_regex.match(line)
        if match:
            factors.append(match.group(1))

    return factors
