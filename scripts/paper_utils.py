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
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    f1_score = 2 * (precision * recall) / (precision + recall)
    return precision, recall, f1_score


def count_on_file(expression, file, compressed=False):
    return int(check_output(f"{'zgrep' if compressed else 'grep'} -c '{expression}' {file}", shell=True))


def parse_links(data):
    vals = data.split(';')
    result = {}
    for val in vals:
        pair = val.split(':')
        try:
            if pair[0] in result.keys():
                result[pair[0]] += f';{pair[1]}'
            else:
                result[pair[0]] = pair[1]
        except:
            pass
    return result


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


def blast_consensus(alignment_file):
    query_to_ref = {}
    ref_to_query = {}
    res = {}
    if not os.path.exists(alignment_file):
        return res
    with open(alignment_file) as file:
        line = file.readline()
        while line:
            line = line.strip('\n')
            line = line.split('\t')
            query_seq = line[0]
            ref_seq = line[1]
            evalue = float(line[-2])
            if query_seq not in query_to_ref:
                if ref_seq not in ref_to_query:
                    query_to_ref[query_seq] = {'ref_seq': ref_seq, 'evalue': evalue}
                    ref_to_query[ref_seq] = {'query_seq': query_seq, 'evalue': evalue}
                else:
                    if ref_to_query[ref_seq]['evalue'] > evalue:
                        ref_to_query[ref_seq] = {'query_seq': query_seq, 'evalue': evalue}
                        query_to_ref[query_seq] = {'ref_seq': ref_seq, 'evalue': evalue}
            else:
                if ref_seq not in ref_to_query:
                    if query_to_ref[query_seq]['evalue'] > evalue:
                        ref_to_query[ref_seq] = {'query_seq': query_seq, 'evalue': evalue}
                        query_to_ref[query_seq] = {'ref_seq': ref_seq, 'evalue': evalue}
                else:
                    if ref_to_query[ref_seq]['evalue'] > evalue:
                        if query_to_ref[query_seq]['evalue'] > evalue:
                            ref_to_query[ref_seq] = {'query_seq': query_seq, 'evalue': evalue}
                            query_to_ref[query_seq] = {'ref_seq': ref_seq, 'evalue': evalue}
            line = file.readline()
    for query_seq in query_to_ref:
        ref_seq = query_to_ref[query_seq]['ref_seq']
        if query_seq == ref_to_query[ref_seq]['query_seq']:
            res[query_seq] = query_to_ref[query_seq]['ref_seq']
    res = pd.DataFrame.from_dict(res, orient='index').reset_index()
    res.columns = ['qseqid', 'Entry']
    res['Entry'] = res['Entry'].apply(lambda x: x.split('|')[1])
    return res


def get_fasta_keys(file):
    return check_output(f"grep '>' {file} | awk '{{print(substr($0, 2))}}'", shell=True).decode('utf8').split('\n')


def parse_gff(filename):
    file = pd.read_csv(filename, sep='\t')
    file.columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    return pd.concat([file, pd.DataFrame(list(file['attribute'].apply(
        lambda x: {part[0]:part[1] for part in [y.split('=') for y in x.split(';')]})))], axis=1)


# is same identification if Jaccard Index / 100 > 0.5
def is_same(list1, list2):
    def jaccard_index(l1, l2):
        return len([ide for ide in l1 if ide in l2]) / len(set(l1 + l2))
    if list1 == [''] or list2 == ['']:
        return False
    return jaccard_index(list1, list2) > 0.5


def get_fasta_ids(filename):
    return [line[1:-1] for line in open(filename) if line.startswith('>')]
