**This contains all required scripts for running the Cell Ranger ARC and the AGGR.
Here you can find a brief description of each script/file:**

**cr_aggr_allsamples.sh**: This script aggregates multiple single-cell multiomic data samples using Cell Ranger ARC from 10x Genomics. It utilizes SLURM for job scheduling, requesting specific computational resources and setting email notifications. The script loads the necessary Cell Ranger ARC module, sets the working directory, and runs the aggregation process with specified parameters, including input libraries, normalization method, reference genome, and resource allocation.

**cr_arc_S180574.sh**: This script processes single-cell multiomic data for a specific sample (ID: S_180574) using Cell Ranger ARC from 10x Genomics. It employs SLURM for job scheduling, specifying resource requirements and setting up email notifications. The script loads the necessary Cell Ranger ARC module, sets the working directory, and runs the counting process with parameters that include the reference genome, input library CSV file, and allocated computational resources.

**cr_arc_S210183.sh**: This script processes single-cell multiomic data for a specific sample (ID: S_210183) using Cell Ranger ARC from 10x Genomics. It employs SLURM for job scheduling, specifying resource requirements and setting up email notifications. The script loads the necessary Cell Ranger ARC module, sets the working directory, and runs the counting process with parameters that include the reference genome, input library CSV file, and allocated computational resources.

**cr_arc_S210414.sh**: This script processes single-cell multiomic data for a specific sample (ID: S_210414) using Cell Ranger ARC from 10x Genomics. It employs SLURM for job scheduling, specifying resource requirements and setting up email notifications. The script loads the necessary Cell Ranger ARC module, sets the working directory, and runs the counting process with parameters that include the reference genome, input library CSV file, and allocated computational resources.

**cr_arc_S211855.sh**: This script processes single-cell multiomic data for a specific sample (ID: S_211855) using Cell Ranger ARC from 10x Genomics. It employs SLURM for job scheduling, specifying resource requirements and setting up email notifications. The script loads the necessary Cell Ranger ARC module, sets the working directory, and runs the counting process with parameters that include the reference genome, input library CSV file, and allocated computational resources.

**libraries_S180574.csv**: This CSV file specifies the libraries used for processing single-cell multiomic data for sample S_180574. It includes paths to FASTQ files, sample identifiers, and the type of library for each dataset. This file will be emplyed for processing of this sample in the cr_arc_S180574.sh script.

**libraries_S210183.csv**: This CSV file specifies the libraries used for processing single-cell multiomic data for sample S_210183. It includes paths to FASTQ files, sample identifiers, and the type of library for each dataset. This file will be emplyed for processing of this sample in the cr_arc_S210183.sh script.

**libraries_S210414.csv**: This CSV file specifies the libraries used for processing single-cell multiomic data for sample S_210414. It includes paths to FASTQ files, sample identifiers, and the type of library for each dataset. This file will be emplyed for processing of this sample in the cr_arc_S210414.sh script.

**libraries_S211855.csv**: This CSV file specifies the libraries used for processing single-cell multiomic data for sample S_211855. It includes paths to FASTQ files, sample identifiers, and the type of library for each dataset. This file will be emplyed for processing of this sample in the cr_arc_S211855.sh script.

**libraries_VHIO.csv**: This CSV file lists the libraries used for processing multiple single-cell multiomic data samples for the VHIO project. It includes paths to FASTQ files, sample identifiers, and the type of library for each dataset. This file is used by the aggregation script (cr_aggr_allsamples.sh) to combine data from multiple samples.
