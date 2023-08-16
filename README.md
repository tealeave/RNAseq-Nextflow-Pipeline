# RNAseq-Nextflow-Pipeline
A Nextflow-based pipeline for processing RNA-seq data, including trimming, alignment with STAR, quantification with RSEM, and differential splicing analysis with rMATS.

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Steps](#pipeline-steps)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## Prerequisites

- Nextflow
- STAR
- RSEM
- rMATS
- trim_galore
- Conda (for rMATS environment)

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/your_username/RNAseq-Nextflow-Pipeline.git
   cd RNAseq-Nextflow-Pipeline
   ```

2. Ensure all required software listed in [Prerequisites](#prerequisites) is installed and available in your `$PATH`.

## Usage

1. Adjust the parameters in the main script to point to your data and reference files.
2. Run the pipeline:
   ```bash
   nextflow run rnaseq_pipe.nf
   ```

## Pipeline Steps

1. **Trimming**: Using `trim_galore` to trim low-quality bases and adapters.
2. **Alignment**: Aligning reads to the reference genome using STAR.
3. **Quantification**: Quantifying gene and isoform expression levels using RSEM.
4. **Differential Splicing Analysis**: (Optional) Using rMATS for differential splicing analysis.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License

[MIT](https://choosealicense.com/licenses/mit/)

## Acknowledgments

- Thanks to the developers of Nextflow, STAR, RSEM, rMATS, and trim_galore for their invaluable tools.
- (Any other acknowledgments or credits you'd like to include)
