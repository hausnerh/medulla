# Utilities for batch processing in medulla projects using jobsub
import sqlite3
from glob import glob
import subprocess
from pathlib import Path

# SQL schema for the configuration table for storing job configurations
SCHEMA_CONFIGURATION = """
CREATE TABLE IF NOT EXISTS configuration (
    jobid INTEGER PRIMARY KEY,
    cfg TEXT NOT NULL
);
"""

# SQL schema for the jobs table for tracking job statuses
SCHEMA_JOBS = """
CREATE TABLE IF NOT EXISTS jobs (
    jobid INTEGER PRIMARY KEY,
    status TEXT,
    FOREIGN KEY (jobid) REFERENCES configuration(jobid)
);
"""

def command(
    curs : sqlite3.Cursor,
    comm : str,
    vals : tuple = None
):
    """
    Execute a command defined in a string using the provided SQLite 
    cursor. Multiple values can be executed if provided as a list.

    Parameters
    ----------
    curs : sqlite3.Cursor
        The SQLite cursor handle.
    comm : str
        The base command.
    vals : tuple
        Values to use as arguments for the sql command (tuple).

    Returns
    -------
    None.
    """
    try:
        if isinstance(vals, list):
            curs.executemany(comm, vals)
        elif vals:
            curs.execute(comm, vals)
        else:
            curs.execute(comm)
    except Exception as e:
        print(e)

def get_samples(
    tml : str,
    batch_size : int
):
    """
    Get the list of samples from the TOML file after filtering the list
    for samples that have been disabled. The batch size is used to
    split samples into multiple separate samples if requested (i.e. for
    processing large samples in smaller chunks).

    Parameters
    ----------
    tml : str
        Path to the TOML file.
    batch_size : int
        Number of files to include in each batch. If <= 0, no batching
        is performed.

    Returns
    -------
    samples : list[dict]
        List of samples that are enabled.
    """
    # Get the initial list of samples from the TOML file that have not
    # been disabled.
    cfg = toml.load(tml)
    samples = cfg.get('sample', [])
    enabled_samples = [s for s in samples if not s.get('disable', False)]

    # Process the samples and batch them if requested.
    batches = []
    for sample in enabled_samples:
        paths = glob(sample['path'])
        if len(paths) == 0:
            raise FileNotFoundError(f"No files found for sample {sample.get('name', '<unknown>')} with path {sample['path']}")
        if batch_size is None or batch_size <= 0:
            batches.append(sample)
        else:
            for i in range(0, len(paths), batch_size):
                batch_paths = paths[i:i+batch_size]
                if len(batch_paths) == 0:
                    continue
                new_sample = sample.copy()
                new_sample['path'] = batch_paths
                batches.append(new_sample)

    # Return the list of enabled samples.
    return batches

def create_new_project(
    project_dir : Path,
    tml : str,
    batch_size : int
):
    """
    Create a new project directory with the necessary subdirectories
    and a SQLite database to manage the project. Each sample in the
    TOML file is added as a separate job in the database, with the
    configuration modified to include only that sample.

    Parameters
    ----------
    project_dir : Path
        Path to the base directory for the job directory.
    tml : str
        Path to the TOML file containing the configuration.
    batch_size : int
        Number of files to process in each batch.

    Returns
    -------
    None.
    """
    # Create the project directory and a subdirectory for job output,
    # if they do not already exist.
    os.makedirs(project_dir, exist_ok=True)
    os.makedirs(project_dir / 'output', exist_ok=True)

    # Connect to the project database. If the database does not exist,
    # it will be created. If the project database does already exist,
    # throw an error because we do not want to overwrite an existing
    # project.
    if (project_dir / 'project.db').exists():
        raise FileExistsError(f"Project database {project_dir / 'project.db'} already exists.")
    conn = sqlite3.connect(project_dir / 'project.db')
    curs = conn.cursor()
    command(curs, SCHEMA_CONFIGURATION)
    command(curs, SCHEMA_JOBS)
    conn.commit()

    # Load the TOML file and get the samples.
    cfg = toml.load(tml)
    samples = get_samples(tml, batch_size)

    # Form a "batch" config for each sample: i.e., each sample gets a
    # copy of the TOML configuration with the [[tree]] list preserved,
    # the [general] section modified to set the 'output' key to its
    # base name plus a batch suffix, and the singular [[sample]]
    # section corresponding to the sample.
    base = cfg['general']['output']
    ins_configurations = []
    ins_jobs = []
    for si, sample in enumerate(samples):
        job_tml = cfg.copy()
        job_tml['general']['output'] = 'output'
        job_tml['sample'] = [sample,]

        ins_configurations.append((si, toml.dumps(job_tml),))
        ins_jobs.append((si, 'pending'))

    # Insert the job configuration into the database.
    command(curs, "INSERT INTO configuration (jobid, cfg) VALUES (?, ?)", ins_configurations)
    command(curs, "INSERT INTO jobs (jobid, status) VALUES (?, ?)", ins_jobs)
    conn.commit()

    conn.close()

def check_project_status(
    project_dir : str,
):
    """
    Check the status of the project by inspecting the job output in the
    project directory.

    Parameters
    ----------
    project_dir : str
        Path to the base directory for the job directory.

    Returns
    -------
    None.
    """
    # Check if the project database exists.
    if not (project_dir / 'project.db').exists():
        raise FileNotFoundError(f"Project database {project_dir / 'project.db'} does not exist.")
    
    # Copy the project database locally to dodge dcache issues.
    subprocess.run(['cp', project_dir / 'project.db', './project.db'], check=True)
    conn = sqlite3.connect('./project.db')
    curs = conn.cursor()

    # Get the list of job outputs in the output directory.
    output_files = glob(str(project_dir / 'output' / 'output_jobid*.root'))
    completed_jobs = [int(Path(f).stem.split('jobid')[-1]) for f in output_files]
    ins = [('completed', jid) for jid in completed_jobs]
    command(curs, "UPDATE jobs SET status = ? WHERE jobid = ?", ins)
    conn.commit()
    conn.close()

    # Replace the project database copy with the updated version.
    subprocess.run(['mv', './project.db', project_dir / 'project.db'], check=True)

    print(f"[INFO] -- Found {len(completed_jobs)} completed jobs.")

def launch_jobsub(
    project_dir : str,
    njobs : int = -1,
):
    """
    Launch jobs using jobsub for the given project directory. If njobs
    is provided, only that many jobs will be launched.

    Parameters
    ----------
    project_dir : str
        Path to the base directory for the job directory.
    njobs : int
        Number of jobs to launch. If None, launch all pending jobs.

    Returns
    -------
    None.
    """
    # Check if the project database exists.
    if not (project_dir / 'project.db').exists():
        raise FileNotFoundError(f"Project database {project_dir / 'project.db'} does not exist.")

    # Copy the project database locally to dodge dcache issues.
    subprocess.run(['cp', project_dir / 'project.db', './project.db'], check=True)
    conn = sqlite3.connect('./project.db')
    curs = conn.cursor()

    # Get the list of pending jobs.
    command(curs, "SELECT jobid FROM jobs WHERE status = 'pending'")
    pending_jobs = [row[0] for row in curs.fetchall()]
    conn.close()

    # Do some checking that the request is sane. Naturally, if there
    # are no pending jobs, there is nothing to launch. Similarly, if
    # the user requested more jobs than are pending, just launch all
    # of the pending jobs.
    if len(pending_jobs) == 0:
        print("[INFO] -- No pending jobs to launch.")
        return
    else:
        print(f"[INFO] -- Found {len(pending_jobs)} pending jobs.")
    if njobs > len(pending_jobs):
        njobs = len(pending_jobs)
        print(f"[INFO] -- Requested number of jobs exceeds pending jobs. Preparing {njobs} jobs instead.")
    if njobs == -1:
        njobs = len(pending_jobs)
        print(f"[INFO] -- No job count specified. Preparing all {njobs} pending jobs.")

    # Form the jobsub command to launch the jobs.
    cmd = [
        'jobsub_submit',
        '-G', 'sbnd',
        '-N', str(njobs),
        '--memory=1800MB',
        '--disk=10GB',
        '--expected-lifetime=1h',
        '--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE',
        "--append_condor_requirements='(TARGET.HAS_Singularity==true)'",
        '--singularity-image=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest',
        f'file://{Path(__file__).resolve().parent / "submit.sh"}',
        '--',
        f'--project-dir={project_dir.resolve()}',
    ]
    print(f"[INFO] -- Launching {njobs} jobs with command: {' '.join(cmd)}")

    # Query the user to confirm that they want to launch the jobs.
    resp = input("Confirm job launch? [Y/N] ")
    if resp.lower() != 'y':
        print("[INFO] -- User aborted job launch.")
        return

    # Launch the jobs. If the command raises an "ExpiredSignatureError"
    # exception, it likely means that the user's token has expired and
    # they need to run `htgettoken` to refresh it. The exception is
    # printed to stdout by jobsub, so we just need to catch it and
    # print a more user-friendly message.
    try:
        out = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        if 'ExpiredSignatureError' in (output := e.stderr.strip()):
            print("[ERROR] -- Job submission failed due to expired token. Please run `htgettoken` to refresh your token and try again.")
        else:
            print(f"[ERROR] -- Job submission failed with error: {output}")
        return
    
    #stdout = out.stdout.strip()
    #print(stdout)
    print(f"[INFO] -- Launched {njobs} jobs.")