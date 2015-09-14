import sys
import os
import re
import math
import utils

import logging
logger = logging.getLogger('Linalg')
cmd_logger = logging.getLogger('cmd')

# Exact configuration here will depends on instance/hardware type.
def run_rankfile(linalg_params):
    logger.info("--- Generating rankfile ---")
    
    machines = linalg_params['machines']
    num_of_mpi = linalg_params['mpi_rows'] * linalg_params['mpi_cols']
    num_of_mach = len(machines)
    num_of_sock = linalg_params['phys_socks_per_machine']
    num_of_cores_per_sock = linalg_params['phys_core_per_sock']
    jobs_assigned_to_mach = 0

    with open(linalg_params['rankfile'], 'wt', encoding='utf-8') as rfile:

        for mach_no in range(0, num_of_mach):
            if mach_no < num_of_mpi % num_of_mach:
                num_of_jobs = num_of_mpi // num_of_mach + 1
            else:
                num_of_jobs = num_of_mpi // num_of_mach

            cores_unassigned = num_of_cores_per_sock * num_of_sock
            socket_counter = {}
            for sock in range(0, num_of_sock):
                socket_counter[sock] = 0

            for job_id in range(0, num_of_jobs):
                rank_no = jobs_assigned_to_mach + job_id
                sock_no = job_id % num_of_sock
                start_core = socket_counter[sock_no]
                cores_to_use = int(math.ceil(cores_unassigned // (num_of_jobs - job_id)))
                end_core = socket_counter[sock_no] + cores_to_use - 1

                # Case for socket splitting
                if end_core >= num_of_cores_per_sock:
                    core_needed = cores_to_use
                    slot_str = ""
                    while core_needed > 0:
                        sock = min(socket_counter, key=socket_counter.get)
                        core_use = (num_of_cores_per_sock - socket_counter[sock] if core_needed >= num_of_cores_per_sock - socket_counter[sock] else core_needed)
                        core_needed -= core_use
                        start_core = socket_counter[sock]
                        end_core = socket_counter[sock] + core_use - 1
                        slot_str += ("{sock}:{start}-{end},"
                            .format(sock=sock, start=socket_counter[sock], end=end_core))
                        socket_counter[sock] += core_use
                    slot_str = slot_str[0:-1]
                    rfile.write("rank {n}={mach} slot={slot}\n"
                        .format(n=rank_no, mach=machines[mach_no], slot=slot_str))
                    cores_unassigned -= cores_to_use
                    continue

                rfile.write("rank {n}={mach} slot={sock}:{start}-{end}\n"
                    .format(n=rank_no, mach=machines[mach_no], sock=sock_no, start=start_core, end=end_core))

                socket_counter[sock_no] += cores_to_use
                cores_unassigned -= cores_to_use
            jobs_assigned_to_mach += num_of_jobs

    logger.info("--- End of generating rankfile ---")

def run_linalg(linalg_params):
    logger.info("--- Beginning MSieve linear algebra ---")

    linalg_cmd = "mpirun -np " + str(linalg_params['mpi_rows'] * linalg_params['mpi_cols'])
    linalg_cmd += " -H " + ",".join(linalg_params['machines'])
    linalg_cmd += " -rf " + linalg_params['rankfile']
    linalg_cmd += " " + os.path.join(linalg_params['msievedir'], 'msieve')
    linalg_cmd += " -nf " + linalg_params['fb_path']
    linalg_cmd += (" -nc2 \"mpi_nrows={rows} mpi_ncols={cols} target_density={td}\""
        .format(rows=linalg_params['mpi_rows'], cols=linalg_params['mpi_cols'], td=linalg_params['target_density']))
    linalg_cmd += " -v -t " + str(linalg_params['threads'])
    linalg_cmd += " -l " + linalg_params['log_path']
    linalg_cmd += " -s " + linalg_params['dat_path']
    linalg_cmd += " " + str(linalg_params['N'])

    cmd_logger.info(linalg_cmd)
    stdout, stderr, ret = utils.run_command(linalg_cmd, include_stdout=True, include_stderr=True, include_returncode=True, logger=logger)

    if ret != 0:
        logger.error("Received error code " + str(ret) + " from Msieve linear algebra. Exiting...")
        sys.exit(1)

    logger.info("--- End of MSieve linear algebra ---")

def run(parameters):

    linalg_paths = ['tasks', 'msieve', 'linalg']
    linalg_keys = {
            "N": int,
            "msievedir": str, 
            "mpi": str, 
            "hosts": str, 
            "target_density": int,
            "phys_socks_per_machine": int, 
            "phys_core_per_sock": int, 
            "threads_per_core": int, 
            "threads": int,
            "rankfile": str, 
            "fb_path": str,
            "log_path": str,
            "dat_path": str,
            }
    linalg_params = parameters.myparams(linalg_keys, linalg_paths)
    linalg_params['machines'] = [ m.strip() for m in linalg_params['hosts'].split(',') if len(m) > 0 ]
    linalg_params['mpi_rows'], linalg_params['mpi_cols'] = [ int(x) for x in linalg_params['mpi'].split("x") ]

    # Create a rankfile based on current mpi configuration
    run_rankfile(linalg_params)

    # Run linear algebra
    run_linalg(linalg_params)
