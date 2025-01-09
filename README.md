# tumorboard WGS DNA and total RNA pipeline


## General info:
Requires a samplesheet containing 5 columns in specific order (tab separated), without headerline:
1) caseID, 2) NPN normal, 3) NPN tumor , 4) NPN tumor RNA, 5) PCGR tumor value


The PCGR tumorcodes can be found at the end of this readme.

Example samplesheet:

    johnDoe 112217976652	111184925465	111184925473    23

The script will automatically look for fastq or cram files in subfolders at the KGVejle dataarchive location. Theres no need to copy or move any input data.

The user can point to a specific folder containing input data using the --fastq or --cram option. 
This is only needed if input data exists outside the data archive (e.g. if data are in personal folders or stored at other KG Vejle servers).

Using the --default option will search for DNA CRAM files (which should be generated automatically by KG Vejle standard preprocessing pipelines), and for RNA-seq fastq files.


## Usage examples:

Analyze the cases in samplesheet.txt, starting with CRAM files (WGS DNA), and fastq files (RNA-seq). let the script find the relevant files at the archive:

        nextflow run KGVejle/tumorboard -r main --samplesheet /path/to/samplesheet.txt --default

Analyze the cases in samplesheet.txt, starting with CRAM files, manually select folder with input data (CRAM):

        nextflow run KGVejle/tumorboard -r main --samplesheet /path/to/samplesheet.txt --cram /path/to/cramfolder/

Analyze the cases in samplesheet.txt, starting with Fastq files, let the script find the relevant fastqfiles at the archive:

        nextflow run KGVejle/tumorboard -r main --samplesheet /path/to/samplesheet.txt --fastqInput

Analyze the cases in samplesheet.txt, starting with Fastq files, manually select folder with input data (Fastq):

        nextflow run KGVejle/tumorboard -r main --samplesheet /path/to/samplesheet.txt --fastq /path/to/fastqfolder/



## Main options:

IMPORT FROM MAIN WHEN DONE!!!

        
## PCGR tumor codes:
    1=Adrenal Gland
    2=Ampulla of Vater
    3=Biliary Tract
    4 = Bladder/Urinary Tract
    5 = Bone
    6 = Breast
    7 = Cervix
    8 = CNS/Brain
    9 = Colon/Rectum
    10 = Esophagus/Stomach
    11 = Eye
    12 = Head and Neck
    13 = Kidney
    14 = Liver
    15 = Lung
    16 = Lymphoid
    17 = Myeloid
    18 = Ovary/Fallopian Tube
    19 = Pancreas
    20 = Peripheral Nervous System
    21 = Peritoneum
    22 = Pleura
    23 = Prostate
    24 = Skin
    25 = Soft Tissue
    26 = Testis
    27 = Thymus
    28 = Thyroid
    29 = Uterus
    30 = Vulva/Vagina

