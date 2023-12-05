![](zz_MolFormatic.png)

Compound database tool, for storing, analyzing and handling bio-assay data and compound information.

### Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)

## Installation

This is a standalone software and have its own GUI. 
Clone the repo, install the requirements and run main.py 


### Example installation steps
git clone https://github.com/ZexiDilling/structure_search.git
cd structure_search
pip install -r requirements.txt

## Usage

The software is build up with 4 main modules. 

1. Compound database and management
   1. A sql database where you can store compound data
   2. In this there are build some chemoinformatic based on RDKIT for doing structure searching and comparison
2. Bio assay data handling
   1. Analyse bio assays
      1. single point 
      2. dose response
   2. The analysing of the data is based on platereader data from a Spark platereader
   3. linking compounds are based on "worklist" for an ECHO-liquid handler
3. LCMS data for checking purity and stability of compounds
   1. The analysing is based on raw data output from a shimadzu LCMS
4. A logistic modul
   1. A sql database that handles basic information


## Examples

missing description and pictures


## Contributing

For now, it is not possible

## license
...