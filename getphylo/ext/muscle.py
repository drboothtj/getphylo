'''
Runs MUSCLE on a provided fasta file.

Functions:
    run_muscle(filename, outname=None) -> None
    get_muscle_version() -> float
'''
import re
import subprocess
import logging
from getphylo.utils import io

def get_muscle_version() -> float:
    '''
    get the muscle version from the command line
        arguments:
            None
        returns
            version_number:
                a float reprisenting the first two parts of the MUSCLE version number
    '''
    with subprocess.Popen(
        ["muscle", "-version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ) as process:
        out, _ = process.communicate()
    # only the first line matters
    version = out.decode().splitlines()[0]
    # the second chunk is all that's relevant
    # e.g. "MUSCLE v3.8.1551" ... or muscle "5.1.linux64 ..."
    version = version.split()[1].lower()
    # remove the leading 'v' if present
    version = version.lstrip("v")
    # grab the first bit to floatify
    try:
        version_number = float(re.search(r"\d\.\d+", version)[0])
        return version_number
    except TypeError as error:
        raise RuntimeError("cannot determine version of MUSCLE") from error


def run_muscle(filename: str, outname=None) -> None:
    '''
    Run MUSCLE aligner on protein fasta file.
        Arguments:
            filename: path to unaligned sequences
            outname: path for the alignment
        Returns:
            None
    '''
    if outname is None:
        outname = "aligned_" + filename
    args = ["muscle"]
    # change the argument format depending on the version of MUSCLE
    # also, MUSCLE 5 is much slower than previous versions so print a warning!
    if get_muscle_version() >= 5.0:
        logging.warning(
            'You are using a MUSCLE version 5 or later. '
            'Be aware that MUSCLE 5 is much slower than previous versions.'
        )
        args.extend(["-align", filename, "-output", outname])
    else:
        args.extend(["-in", filename, "-out", outname])
    command = " ".join(args)
    io.run_in_command_line(command)
