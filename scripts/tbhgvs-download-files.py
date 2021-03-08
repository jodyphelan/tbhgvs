import sys
import argparse
import tbhgvs
import os


def main(args):
    tbhgvs.download_files(directory=args.dir)
    sys.stdout.write("Files downloaded successfully\n")

# Set up the parser
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--dir',default="%s/.tbhgvs/" % os.path.expanduser("~"),type=str,help='Directory to store files')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
