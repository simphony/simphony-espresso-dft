""" LAMMPS Process
This module provides a way to run the lammps or liggghts process
"""

import os
import subprocess
import logging

from simespresso.io_QE.data_handler import  qe_data_handler

class QeProcess(object):
    """ Class runs the lammps/liggghts program
    Parameters
    ----------
    lammps_name : str
        name of LAMMPS executable
    log_directory : str, optional
        name of directory of log file ('log.lammps') for lammps.
        If not given, then pwd is where 'log.lammps' will be written.
    Raises
    ------
    RuntimeError
        if Lammps did not run correctly
    """
    def __init__(self, qe_executable="pw.x", log_directory=None):
        self._qe_executable = qe_executable
        self._returncode = 0
        self._stderr = ""
        self._stdout = ""
        self._datahandler = qe_data_handler()

        if log_directory:
            self._log = os.path.join(log_directory, 'log.qe')
        else:
            self._log = 'log.qe'

        # see if lammps can be started
        try:
            self.run(" ")
        except Exception:
            msg = "quantum espresso could not be started."
            if self._returncode == 127:
                msg += " executable '{}' was not found.".format(qe_executable)
            else:
                msg += " stdout/err: " + self._stdout + " " + self._stderr
            raise RuntimeError(msg)

    def run(self,
                input_data_file=input_data_filename,
                output_data_file=output_data_filename,
                BC=_combine(self.BC, self.BC_extension),
                CM=_combine(self.CM, self.CM_extension),
                SP=_combine(self.SP, self.SP_extension)):
        """Run engine with a set of commands
        Parameters
        ----------
        commands : str
            set of commands to run
        Raises
        ------
        RuntimeError
            if Lammps did not run correctly
        """
        logging.debug('starting qe engine')

        self._datahandler.write_espresso_input_file(input_data_filename
        logging.debug('path to espresso:'+self._qe_executable)
        if not which(self._qe_executable):
            logging.debug('no path to espresso')
            raise ValueError(
                'espresso command not found (looking for '
                            + self.path_to_espresso+')')

        if self._datahandler.mpi:
            command = 'mpirun -np ' + str(self.mpi_Nprocessors) + ' ' + \
                     path + ' < ' + pwname + ' > ' \
                      + self.datahandler.output_filename
        else:
            command = self.path_to_espresso + ' < ' + pwname + ' > ' \
                      + self.datahandler.output_filename
        logging.debug('attempting to run command: ' + command)
        try:
            subprocess.check_call(command, shell=True,
                stdout=subprocess.PIPE).stdout.read()
#            subprocess.Popen(command, shell=True,
#                             stdout=subprocess.PIPE).stdout.read()
        except:
            e = sys.exc_info()[0]
            logging.debug('espresso command gave error %s' % e)
            raise ValueError('espresso command gave error %s' % e)


'''
        proc = subprocess.Popen(
            [self._qe_executable, '-log', self._log], stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self._stdout, self._stderr = proc.communicate(commands)
        self._returncode = proc.returncode

        if self._returncode != 0 or self._stderr:
            msg = "quantum espresso ('{}') did not run correctly. ".format(
                self._qe_executable)
            msg += "Error code: {} ".format(proc.returncode)
            if self._stderr:
                msg += "stderr: \'{}\n\' ".format(self._stderr)
            if self._stdout:
                msg += "stdout: \'{}\n\'".format(self._stdout)
            raise RuntimeError(msg)
'''


def which(name):
    found = 0
    for path in os.getenv("PATH").split(os.path.pathsep):
        full_path = path + os.sep + name
        if os.path.exists(full_path):
            """
            if os.stat(full_path).st_mode & stat.S_IXUSR:
                found = 1
                print(full_path)
            """
            found = 1
            print(full_path)
    # Return a UNIX-style exit code so it can be checked by calling scripts.
    # Programming shortcut to toggle the value of found: 1 => 0, 0 => 1.
    sys.exit(1 - found)