import sys
import argparse
import tbhgvs

def hgvs(args):
    ref_db = tbhgvs.reference_db()
    mutations = ref_db.hgvs2genome(args.mutation,args.gene)
    print("Pos\tRef\tAlt")
    for mut in mutations:
        positions = ",".join([str(x["pos"]) for x in mut])
        refs = ",".join([x["ref"] for x in mut])
        alts = ",".join([x["alt"] for x in mut])
        print("%s\t%s\t%s" % (positions,refs,alts))

def genome(args):
    ref_db = tbhgvs.reference_db()
    mutations = []
    positions = args.pos.split(",")
    refs = args.ref.split(",")
    alts = args.alt.split(",")
    mutations = [{"pos":x[0],"ref":x[1],"alt":x[2]} for x in zip(positions,refs,alts)]

    mut = ref_db.genome2hgvs(mutations)
    for key,val in mut.items():
        print("%s: %s" % (key,val))


# Set up the parser
parser = argparse.ArgumentParser(description='tbprofiler script',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_sub = subparsers.add_parser('hgvs', help='Convert from hgvs to genome', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--gene',help='Mutation panel name (e.g. rpoB)',required=True)
parser_sub.add_argument('--mutation',help='Mutation (e.g. p.Ser450Leu)',required=True)
parser_sub.set_defaults(func=hgvs)

parser_sub = subparsers.add_parser('genome', help='Convert from hgvs to genome', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--pos',help='Positions',required=True)
parser_sub.add_argument('--ref',help='Reference alleles',required=True)
parser_sub.add_argument('--alt',help='Alternate alleles',required=True)
parser_sub.set_defaults(func=genome)

args = parser.parse_args()
if vars(args)=={}:
    parser.print_help(sys.stderr)
else:
    args.func(args)
