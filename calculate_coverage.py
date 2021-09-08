# samtools index bam for IGV
import subprocess
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--path', type=str, help='path to sorted bam file')
parser.add_argument('-o', '--output', type=str, help='output file')
args = parser.parse_args()

bamPath = args.path
covFile = args.output

fileList = {}

for f in os.listdir(bamPath):
    if f.endswith('.bam'):
        filePath = os.path.join(bamPath, f)
        sample = os.path.splitext(f)[0]
        fileList[sample] = filePath

bamListFile = os.path.join(bamPath, 'bamlist.txt')
with open(bamListFile, 'w') as listFile:
    for sample in fileList:
        listFile.write(fileList[sample] + '\n')

args = ['samtools', 'depth',
        '-aa',
        '-f', bamListFile
        ]
print('Running...\n', ' '.join(args))
p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
with open(covFile, 'w') as outputHandle:
    # write header
    outputHandle.write(
        'chr\tpos\t' + '\t'.join([sample for sample in fileList]) + '\n')
    # write output
    for line in p.stdout:
        data = line.decode('utf-8')
#        print(data, end='')
        outputHandle.write(data)
