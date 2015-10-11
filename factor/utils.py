import os
import shlex
import subprocess
import time
import logging
from cadotask import Polynomials
import cadoprograms
import cadoparams
import datetime
import fcntl
import pickle
import copy

import threading

def str_time(t):
    m, s = divmod(t, 60)
    h, m = divmod(m, 60)
    return "%02d:%02d:%02d" % (h, m, s)
    # return time.strftime('%H:%M:%S', time.mktime(t))

def parse_time(s):
    h, m, s = s.split(':')
    return 60 * 60 * int(h) + 60 * int(m) + int(s)

checkpoint_file = None
checkpoint_lock = threading.Lock()

def init_checkpoint(resume_file):
    global checkpoint_file
    checkpoint_file = resume_file
    if not os.path.isfile(checkpoint_file):
        set_checkpoint({})
    checkpoint = get_checkpoint()
    if checkpoint.get('stage', None) == None:
        set_checkpoint({'stage': 0})

def clear_checkpoint():
    set_checkpoint({})

def update_checkpoint(update):
    checkpoint = get_checkpoint()
    checkpoint.update(update)
    set_checkpoint(checkpoint)

def set_checkpoint(checkpoint):
    with checkpoint_lock:
        with open(checkpoint_file, 'wb') as f:
            pickle.dump(checkpoint, f)

def get_checkpoint():
    with checkpoint_lock:
        with open(checkpoint_file, 'rb') as f:
            return pickle.load(f)

def non_block_read(output):
    fd = output.fileno()
    fl = fcntl.fcntl(fd, fcntl.F_GETFL)
    fcntl.fcntl(fd, fcntl.F_SETFL, fl | os.O_NONBLOCK)
    try:
        return output.readline()
    except:
        return b""


slurm_logfile = '/var/log/slurm/accounting'
def slurm_cpu_time_start():
    offset = os.path.getsize(slurm_logfile)
    return offset

def slurm_cpu_time_end(offset, outputfile):
    jobids = set()
    with open(slurm_logfile) as f:
        f.seek(offset)
        for line in f:
            # Parse the job id from the line
            jobid = line.split()[0].split('.')[0]
            jobids.add(jobid)
    cputime_total = 0
    batch = set()
    count = 0
    total = len(jobids)
    while len(jobids) > 0:
        batch.add(jobids.pop())
        count += 1
        # max subprocess buffer length is 32,768 characters, so run this in batches
        if len(batch) >= 1000 or len(jobids) == 0:
            cmd = 'sacct --noheader --parsable2 --format=JobName,JobId,Elapsed,CPUTime -s COMPLETED -j %s -X' % ','.join(batch)
            output = run_command(cmd)
            with open(outputfile, 'a') as f:
                f.write(output)
            for line in output.split('\n'):
                if len(line) > 0: 
                    jobname,jobid,elapsed,cputime = line.split('|')
                    cputime_total += parse_time(cputime)
            batch.clear()
    return cputime_total

def run_command(command, include_stdout=True, include_stderr=False, include_time=False, include_returncode=False, logger=None):
    command_args = shlex.split(command)
    stime = os.times()

    child = None
    stdout = None
    stderr = None
    if logger != None:
        child = subprocess.Popen(command_args, bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout_done = False
        stderr_done = False
        stdout = b""
        stderr = b""
        while not stdout_done or not stderr_done:
            stdout_line = non_block_read(child.stdout)
            if not stdout_line:
                pass 
            elif stdout_line == "":
                stdout_done = True
            else:
                logger.debug("stdout: " + stdout_line.decode('utf-8').strip())
                stdout += stdout_line

            stderr_line = non_block_read(child.stderr)
            if not stderr_line:
                pass
            elif stderr_line == "":
                stderr_done = True
            else:
                logger.debug("stderr: " + stderr_line.decode('utf-8').strip())
                stderr += stderr_line

            if child.poll() != None:
                break
        # Read the remaining output in the pipes
        stdout_remaining = child.stdout.read()
        if stdout_remaining:
            logger.debug("stdout: " + stdout_remaining.decode('utf-8').strip())
            stdout += stdout_remaining

        stderr_remaining = child.stderr.read()
        if stderr_remaining:
            logger.debug("stderr: " + stderr_remaining.decode('utf-8').strip())
            stderr += stderr_remaining

        child.stdout.close()
        child.stderr.close()

    else:
        child = subprocess.Popen(command_args, bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
        (stdout, stderr) = child.communicate()

    ftime = os.times()

    to_return = []
    if include_stdout:
        to_return.append(stdout.decode('utf-8'))
    if include_stderr:
        to_return.append(stderr.decode('utf-8'))
    if include_time:
        to_return.append(ftime[2] + ftime[3] - stime[2] + stime[3])
    if include_returncode:
        to_return.append(child.returncode)
    if len(to_return) == 1:
        return to_return[0]
    else:
        return tuple(to_return)

def get_params(paramfile, options):
    parameters = cadoparams.Parameters()
    parameters.readfile(paramfile)
    parameters.readparams(options)
    return parameters


def write_list_to_file(files_list, root_dir):
    files_name = os.path.join(root_dir, 'filelist.%d' % (time.time()))
    with open(files_name, 'w') as f:
        f.write('\n'.join(files_list) + '\n')
    return files_name

# input is in format outputted by polysel1 and polysel2
def parse_poly(input):
    block = []
    for line in input.split("\n"):
        line = line.strip()
        if line:
            block.append(line)
        else:
            if block:
                try:
                    yield Polynomials(block)
                except Exception as e:
                    pass
            block = []
    if block:
        try:
            yield Polynomials(block)
        except Exception as e:
            pass
