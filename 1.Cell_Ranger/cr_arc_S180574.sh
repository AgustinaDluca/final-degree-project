#!/bin/bash
#SBATCH --job-name=180574cr_arc
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=120gb                    # Job memory request
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --mail-type=FAIL,END
#SBATCH --chdir="/ijc/LABS/MEREU/DATA/aitor/VHIO_multiome/cellranger_arc"
#SBATCH --mail-user=adeluca@carrerasresearch.org



module load cellranger-arc/2.0.2

cellranger-arc count --id=S_180574 \
                       --reference=/ijc/LABS/MEREU/RAW/Cellranger_references/ARC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                       --libraries=/ijc/LABS/MEREU/DATA/aitor/VHIO_multiome/cellranger_arc/libraries_S180574.csv \
                       --localcores=16 \
                       --localmem=120
