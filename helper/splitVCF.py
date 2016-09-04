import argparse,tempfile
from os.path import join,exists
from os import remove,chdir,system,makedirs


def parse_args():
    parser = argparse.ArgumentParser(description="Convert sequence and target for Caffe")

    # Positional (unnamed) arguments:
    parser.add_argument("vcffile",  type=str, help="VCF file to split")
    parser.add_argument("outdir",type=str,help="Output dir")
    return parser.parse_args()

if __name__ == "__main__":

    args = parse_args()

    if exists(args.outdir):
        system('rm -r ' + args.outdir)

    makedirs(args.outdir)


    tmp = tempfile.NamedTemporaryFile().name
    with open(tmp,'w') as fout,open(args.vcffile,'r') as fin:
        for x in fin:
            if x[0] != '#':
                fout.write('%s' % x)

    chdir(args.outdir)
    cmd = ' '.join(['awk -F\"\\t\" \'{ print > (\"chr\" $1 \".vcf\") }\'',tmp])
    system(cmd)

    remove(tmp)
