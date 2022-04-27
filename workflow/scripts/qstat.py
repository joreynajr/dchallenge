import subprocess as sp
cmd = '/usr/local/bin/qstat'
result = sp.check_output(cmd, shell=True)
print(result)
