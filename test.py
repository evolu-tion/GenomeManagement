import os
import sys

from argparse import ArgumentParser

args = parse_cmdline(sys.argv)

# Process command-line arguments
def parse_cmdline(args):
    """Parse command-line arguments for script."""
    parser = ArgumentParser(prog="average_nucleotide_identity.py")
    parser.add_argument("-o", "--outdir", dest="outdirname",
                        action="store", default=None,
                        help="Output directory")
    parser.add_argument("-i", "--indir", dest="indirname",
                        action="store", default=None,
                        help="Input directory name")
    parser.add_argument("-v", "--verbose", dest="verbose",
                        action="store_true", default=False,
                        help="Give verbose output")
    parser.add_argument("-f", "--force", dest="force",
                        action="store_true", default=False,
                        help="Force file overwriting")
    parser.add_argument("-s", "--fragsize", const="fragsize",
                        action="store_const", default=pyani_config.FRAGSIZE,
                        help="Sequence fragment size for ANIb")
    parser.add_argument("-l", "--logfile", dest="logfile",
                        action="store", default=None,
                        help="Logfile location")
    parser.add_argument("--skip_nucmer", dest="skip_nucmer",
                        action="store_true", default=False,
                        help="Skip NUCmer runs, for testing " +
                        "(e.g. if output already present)")
    parser.add_argument("--skip_blastn", dest="skip_blastn",
                        action="store_true", default=False,
                        help="Skip BLASTN runs, for testing " +
                        "(e.g. if output already present)")
    parser.add_argument("--noclobber", dest="noclobber",
                        action="store_true", default=False,
                        help="Don't nuke existing files")
    parser.add_argument("-g", "--graphics", dest="graphics",
                        action="store_true", default=False,
                        help="Generate heatmap of ANI")
    parser.add_argument("--gformat", dest="gformat",
                        action="store", default="pdf,png,eps",
                        help="Graphics output format(s) [pdf|png|jpg|svg]")
    parser.add_argument("--gmethod", dest="gmethod",
                        action="store", default="mpl",
                        help="Graphics output method [mpl|R]")
    parser.add_argument("--labels", dest="labels",
                        action="store", default=None,
                        help="Path to file containing sequence labels")
    parser.add_argument("--classes", dest="classes",
                        action="store", default=None,
                        help="Path to file containing sequence classes")
    parser.add_argument("-m", "--method", dest="method",
                        action="store", default="ANIm",
                        help="ANI method [ANIm|ANIb|ANIblastall|TETRA]")
    parser.add_argument("--scheduler", dest="scheduler",
                        action="store", default="multiprocessing",
                        help="Job scheduler [multiprocessing|SGE]")
    parser.add_argument("--maxmatch", dest="maxmatch",
                        action="store_true", default=False,
                        help="Override MUMmer to allow all NUCmer matches")
    parser.add_argument("--nucmer_exe", dest="nucmer_exe",
                        action="store", default=pyani_config.NUCMER_DEFAULT,
                        help="Path to NUCmer executable")
    parser.add_argument("--blastn_exe", dest="blastn_exe",
                        action="store", default=pyani_config.BLASTN_DEFAULT,
                        help="Path to BLASTN+ executable")
    parser.add_argument("--makeblastdb_exe", dest="makeblastdb_exe",
                        action="store",
                        default=pyani_config.MAKEBLASTDB_DEFAULT,
                        help="Path to BLAST+ makeblastdb executable")
    parser.add_argument("--blastall_exe", dest="blastall_exe",
                        action="store", default=pyani_config.BLASTALL_DEFAULT,
                        help="Path to BLASTALL executable")
    parser.add_argument("--formatdb_exe", dest="formatdb_exe",
                        action="store",
                        default=pyani_config.FORMATDB_DEFAULT,
                        help="Path to BLAST formatdb executable")
    parser.add_argument("--write_excel", dest="write_excel",
                        action="store_true",
                        default=False,
                        help="Write Excel format output tables")
    return parser.parse_args()
