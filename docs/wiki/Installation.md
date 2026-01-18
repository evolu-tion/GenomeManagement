# Installation

GenomeManagement requires Python 3.7 or higher.

## Install from PyPI

Once released, you can install the package directly from PyPI:

```bash
pip install genomemanagement
```

## Install from GitHub

The easiest way to install the package is directly from the GitHub repository:

```bash
pip install git+https://github.com/evolu-tion/GenomeManagement.git
```

## Install from Source

You can clone the repository and install it locally.

### Standard Installation

Use this for normal usage:

1.  Clone the repository:
    ```bash
    git clone https://github.com/evolu-tion/GenomeManagement.git
    ```
2.  Navigate to the directory:
    ```bash
    cd GenomeManagement
    ```
3.  Install the package:
    ```bash
    pip install .
    ```

### Development Installation

Use this if you plan to modify the code:

```bash
pip install -e .
```

## Verify Installation

After installation, verify that the package is correctly installed by checking the version of one of the tools:

```bash
promoter-retrieve --version
```

You should see output similar to: `GenomeManagement_v2.0.0`
