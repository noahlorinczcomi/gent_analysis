import subprocess
import os
import argparse
# python is only used to create the slurm file
os.chdir('/home/lorincn/isilon/Cheng-Noah/GenT_results/gent_specific_traits_ondemand')

parser=argparse.ArgumentParser(prog='GenT',description='Perform genome-wide GenT analysis.')
##########
# defaults
ld_default='/home/lorincn/beegfs/lorincn/data/reference_panels/1kg.v3/EUR'
mem_gb=50
##########
parser.add_argument('--gwas',action='store',type=str,help='(Required) filepath with .gz extension to phentoype GWAS.')
parser.add_argument('--jobname',action='store',type=str,default='gent_on_demand',help='SLURM job name.')
parser.add_argument('--ldref',action='store',type=str,default=ld_default,help='PLINK-formatted LD reference (without extension).')
parser.add_argument('--memgb',action='store',type=int,default=mem_gb,help='Gb or memory to request (50gb by default).')
args=parser.parse_args()
# edit .R file that performs GenT
## edits lines 5 (GWAS filepath) and 6 (LD reference)
### read file
with open('gent_specific_trait.R','r') as f:
    data=f.readlines()

### edit file
data[4]='gwas_fp='+'"'+args.gwas+'"\n'
data[5]='ldref='+'"'+args.ldref+'"\n'
### write file
with open('temp.R','w') as f:
    f.writelines(data)

# create .slurm file that runs the job
with open('temp.slurm','w') as f:
    f.write('''#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lorincn@ccf.org
#SBATCH --job-name='''+str(args.jobname)+'''
#SBATCH --mem='''+str(args.memgb)+'''gb
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --partition=defq

cd /home/lorincn/isilon/Cheng-Noah/GenT_results/gent_specific_traits_ondemand
module load R/4.3.3 plink/1.90
Rscript temp.R
''')

out=subprocess.call(['sbatch', 'temp.slurm'],stdout=subprocess.DEVNULL,stderr=subprocess.STDOUT)
