'''
Make DIAMOND databases and run searches using DIAMOND

Functions:
make_diamond_database(filename)
run_diamond_search(filename)

'''
from getphylo import console

def make_diamond_database(filename: str, dmnd_database=None):
    '''Create a DIAMOND database from a fasta file'''
    if dmnd_database is None:
        database = " --db " + filename.split('.')[0] + ".dmnd"
    else:
        database = " --db " + dmnd_database
    makedb = "diamond makedb"
    input_ = " --in " + filename + ""
    command = makedb + database + input_
    console.run_in_command_line(command)

def run_diamond_search(filename: str, dmnd_database=None, outname=None):
    '''Run BLASTP through DIAMOND'''
    if dmnd_database is None:
        database = " --db " + filename.split('.')[0] + ".dmnd"
    else:
        database = " --db " + dmnd_database
    if outname is None:
        output = " --out " + dmnd_database.split('.')[0] + ".tsv"
    else:
        output = " --out " + outname
    blastp = "diamond blastp"
    query = " --query " + filename
    outfmt = " --outfmt 6 qseqid sseqid pident"
    command = blastp + database + query + output + outfmt
    console.run_in_command_line(command)
