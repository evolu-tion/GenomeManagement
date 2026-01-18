# Publishing on Bioconda

Bioconda is a channel for the conda package manager specializing in bioinformatics software. Publishing there makes your tool easily installable via `conda install -c bioconda genomemanagement`.

Unlike PyPI, you do **not** push directly to Bioconda. Instead, you submit a "Recipe" (a configuration file) via a Pull Request to the [bioconda-recipes](https://github.com/bioconda/bioconda-recipes) repository.

## Prerequisites

1.  **Release on PyPI**: Ensure your package (v2.0.0) is successfully published on PyPI.
2.  **GitHub Account**: You need a GitHub account to submit a Pull Request.

## Step-by-Step Guide

### 1. Get the SHA256 Checksum

Once your package is on PyPI, find the SHA256 checksum for the source tarball (`.tar.gz`).
You can find this on the PyPI "Download files" page for your release, or run:

```bash
curl -sL https://pypi.io/packages/source/g/genomemanagement/genomemanagement-2.0.0.tar.gz | shasum -a 256
```

### 2. Fork bioconda-recipes

Go to [bioconda-recipes](https://github.com/bioconda/bioconda-recipes) and click **Fork** to create a copy in your account.

### 3. Create the Recipe

In your forked repository, create a new directory: `recipes/genomemanagement`.
Inside that directory, create a file named `meta.yaml` with the following content:

```yaml
{% set name = "genomemanagement" %}
{% set version = "2.0.0" %}

package:
  name: {{ name|lower }}
  version: {{ version }}

source:
  url: https://pypi.io/packages/source/{{ name[0] }}/{{ name }}/{{ name }}-{{ version }}.tar.gz
  sha256: <PASTE_YOUR_SHA256_HERE>

build:
  number: 0
  noarch: python
  script: {{ PYTHON }} -m pip install . -vv
  run_exports:
    - {{ pin_subpackage('genomemanagement', max_pin="x") }}
  entry_points:
    - promoter-retrieve = genomemanagement.promoter_retrieve:main
    - proteins-retrieve = genomemanagement.proteins_retrieve:main
    - get-genome-statistic = genomemanagement.get_genome_statistic:main
    - plantpan2 = genomemanagement.PlantPAN2:main
    - plantpan3 = genomemanagement.PlantPAN3:main
    - meme-xml2bed = genomemanagement.memeXml2bed:main

requirements:
  host:
    - python >=3.7
    - pip
    - setuptools
  run:
    - python >=3.7

test:
  imports:
    - genomemanagement
  commands:
    - promoter-retrieve --help
    - proteins-retrieve --help
    - get-genome-statistic --help

about:
  home: https://github.com/evolu-tion/GenomeManagement
  summary: "Genome Management Package for bioinformatics analysis"
  license: MIT
  license_file: LICENSE

extra:
  recipe-maintainers:
    - evolu-tion  # Your GitHub username
```

### 4. Submit Pull Request

1.  Commit the `meta.yaml` file to your fork.
2.  Open a Pull Request against the original `bioconda/bioconda-recipes` repository.
3.  The Bioconda bots will automatically lint and test your recipe.
4.  Once the tests pass, a Bioconda maintainer will merge it.

### 5. Install

After merging, wait a few hours, then your package will be available:

```bash
conda install -c bioconda genomemanagement
```
