'''
Use strings to interact with the terminal.

Functions:
print_to_system(string_to_print: str)
run_in_command_line(command) -> process

'''

from datetime import datetime
import subprocess

def print_to_system(string_to_print: str):
    '''Print a string to the terminal with the current time'''
    now = datetime.now()
    current_time = now.strftime("[%H:%M:%S]: ")
    print(current_time + string_to_print)

def run_in_command_line(command: str):
    '''Convert a string into a command and run in the terminal'''
    command = command.split(" ")
    process = subprocess.Popen(command, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    process.communicate()
    return process
