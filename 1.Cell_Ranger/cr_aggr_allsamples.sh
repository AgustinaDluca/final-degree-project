#!/bin/sh
#SBATCH --job-name=agus_agg_arc
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=120gb                    # Job memory request
#SBATCH --time=72:00:00               # Time limit hrs:min:sec
#SBATCH --mail-type=FAIL,END
#SBATCH --chdir="/ijc/LABS/MEREU/DATA/aitor/VHIO_multiome/cellranger_arc"
#SBATCH --mail-user=adeluca@carrerasresearch.org



module load cellranger-arc/2.0.2

cellranger-arc aggr --id=agus_agg_arc \
                    --csv=/ijc/LABS/MEREU/DATA/aitor/VHIO_multiome/cellranger_arc/libraries.csv \
                    --normalize=depth \
                    --reference=/ijc/LABS/MEREU/RAW/Cellranger_references/ARC/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                    --localcores=16 \
                    --localmem=120

