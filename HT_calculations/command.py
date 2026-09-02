class Command(object):
    """
    copy from the codes in 2DMSID by xlzhang

    Helper class -  run subprocess commands in a different thread with TIMEOUT option.
    From https://gist.github.com/kirpit/1306188
    Based on jcollado's solution:
    http://stackoverflow.com/questions/1191374/subprocess-with-timeout/4825933#4825933
    """
    command = None
    process = None
    status = None
    output, error = '', ''

    def __init__(self, command):
        """
        initialize the object.

        Args:
            command: command to run
        """
        import shlex
        if isinstance(command, str):
            command = shlex.split(command)
        self.command = command

    def run(self, timeout=None, **kwargs):
        """
        Run the command.

        Args:
            timeout (float): timeout
            kwargs (dict)

        Returns:
            (status, output, error)
        """
        import subprocess, threading
        def target(**kwargs):
            try:
                self.process = subprocess.Popen(self.command, **kwargs)
                self.output, self.error = self.process.communicate()
                self.status = self.process.returncode

                # Python3 - need to convert string to bytes
                if isinstance(self.output, bytes):
                    self.output = self.output.decode("utf-8")
                if isinstance(self.error, bytes):
                    self.error = self.error.decode("utf-8")
            except Exception:
                self.status = -1

        # default stdout and stderr
        if 'stdout' not in kwargs:
            kwargs['stdout'] = subprocess.PIPE
        if 'stderr' not in kwargs:
            kwargs['stderr'] = subprocess.PIPE
        # thread
        thread = threading.Thread(target=target, kwargs=kwargs)
        thread.start()
        thread.join(timeout)
        if thread.is_alive():
            self.process.terminate()
            thread.join()
        return self.status, self.output, self.error

    pass
