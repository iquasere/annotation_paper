from subprocess import run, check_output, Popen, PIPE
import os
import pandas as pd


def get_fasta_ids(file):
    return [ide[1:] for ide in check_output(f"grep '>' {file}", shell=True).decode('utf8').split()]


def run_command(bash_command, print_message=True):
    if print_message:
        print(bash_command)
    run(bash_command.split(), check=True)


def run_pipe_command(bashCommand, output='', mode='w', print_message=True):
    if print_message:
        print(bashCommand)
    if output == '':
        Popen(bashCommand, stdin=PIPE, shell=True).communicate()
    elif output == 'PIPE':
        return Popen(bashCommand, stdin=PIPE, shell=True, stdout=PIPE).communicate()[0].decode('utf8')
    else:
        with open(output, mode) as output_file:
            Popen(bashCommand, stdin=PIPE, shell=True, stdout=output_file).communicate()



def parse_blast(file):
    if os.stat(file).st_size != 0:
        blast = pd.read_csv(file, sep='\t', header=None)
        blast.columns = [
            'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
            'bitscore']
        return blast
    return pd.DataFrame(columns=[
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
        'bitscore'])


def calculate_quality_metrics(tp, fp, fn):
    precision = tp/(tp+fp)
    recall = tp/(tp+fn)
    f1_score = 2 * (precision * recall) / (precision + recall)
    return precision, recall, f1_score


def count_on_file(expression, file, compressed=False):
    return int(check_output(f"{'zgrep' if compressed else 'grep'} -c '{expression}' {file}", shell=True))


def parse_mantis_consensus(filename):
    file = open(filename).readlines()
    result = []
    for line in file[1:]:
        line = line.rstrip('\n').split('\t')
        result.append(line[:5] + [';'.join(line[6:])])
    df = pd.DataFrame(result, columns=['Query', 'Ref_Files', 'Ref_Hits', 'Consensus_hits', 'Total_hits', 'Links'])
    df = pd.concat([df, pd.DataFrame.from_records(df['Links'].apply(parse_links))], axis=1)
    del df['Links']
    return df