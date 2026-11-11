import sys
import re

def get_number(string):
    print(int(float(re.findall(r'\d+', string)[0])))
    return int(float(re.findall(r'\d+', string)[0]))

print(len(sys.argv))
print(sys.argv[1:])
print(sorted(sys.argv[1:]))
print(sorted(sys.argv[1:], key=get_number))
print(sys.argv[1:].sort(key=get_number))
