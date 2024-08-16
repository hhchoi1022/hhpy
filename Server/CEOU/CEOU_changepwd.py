import string
import random
import os
from astropy.table import Table
import time
 
new_pw_len = 10 # 새 비밀번호 길이
pw_candidate = string.ascii_letters + string.digits + '!@#$%^'
f = open('/etc/passwd')
lines = f.readlines()
IDlist = []
pwdlist = []
for line in lines:
    new_pw = ""
    keys = line.split(':')
    user_key = keys[-1]
    if user_key.endswith('bash\n'):
        for i in range(new_pw_len):
            new_pw += random.choice(pw_candidate)
        ID = keys[0]
        command = "echo "+new_pw+' | passwd --stdin '+ID
        print(command)
        os.system(command)
        time.sleep(0.5)
        IDlist.append(ID)
        pwdlist.append(new_pw)
        command = ""

result = Table()
result['ID'] = IDlist
result['pwd'] = pwdlist
result.write('temporary.txt',format = 'ascii.fixed_width')

