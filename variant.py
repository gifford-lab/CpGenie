from os.path import join,basename,dirname,exists,realpath
from os import makedirs,system,remove
from tempfile import mktemp,mkdtemp
import numpy as np,multiprocessing as mp,sys,time

cwd = dirname(realpath(__file__))
def getseq(seq,revseq,startpos,endpos,f,wind,anno):
    ### startpos and end pos should be all 1-based
    tchr = '1'
    offset = -2*wind
    cnt = 0
    all_idx = []
    for idx in range(startpos-1,endpos):
        if seq[idx] == 'C' and seq[idx+1] == 'G' and idx != (endpos-1):
            f.write('\t'.join(['>chr'+tchr+':'+str(startpos)+':'+str(endpos)+':'+str(len(seq))+':'+str(idx-wind+offset)+':'+str(idx+wind+offset)+anno,''.join(seq[(idx-wind):(idx+wind+1)])]))
            f.write('\n')
            cnt +=1
            all_idx.append(idx)
        if seq[idx] == 'G' and seq[idx-1] == 'C' and idx != (startpos-1):
            f.write('\t'.join(['>rev:chr'+tchr+':'+str(startpos)+':'+str(endpos)+':'+str(len(seq))+':'+str(idx-wind+offset)+':'+str(idx+wind+offset)+anno,''.join(revseq[(idx-wind):(idx+wind+1)][::-1])]))
            f.write('\n')
            cnt += 1
            all_idx.append(idx)
    return cnt

def slave(args):
    data,wind,index,outdir,anno,cpg_range = args[:]
    ref_seq,alt_seq = data[:]
    dict={'A':'T','C':'G','G':'C','T':'A','N':'N'}
    ref_seq_rev = [dict[x]for x in ref_seq]
    alt_seq_rev = [dict[x]for x in alt_seq]
    ref_2evlfa = join(outdir,'test.refevl.'+str(index))
    alt_2evlfa = join(outdir,'test.altevl.'+str(index))

    ref_len = 0
    alt_len = 0
    ref_seq_len = len(ref_seq)
    alt_seq_len = len(alt_seq)
    with open(ref_2evlfa,'w') as fref,open(alt_2evlfa,'w') as falt:
        ref_len += getseq(ref_seq,ref_seq_rev,wind+1,wind+cpg_range,fref,wind,anno[0])
        alt_len += getseq(alt_seq,alt_seq_rev,wind+1,wind+cpg_range,falt,wind,anno[1])
        ref_len += getseq(ref_seq,ref_seq_rev,ref_seq_len-wind-cpg_range+1,ref_seq_len-wind,fref,wind,anno[0])
        alt_len += getseq(alt_seq,alt_seq_rev,alt_seq_len-wind-cpg_range+1,alt_seq_len-wind,falt,wind,anno[1])
    try:
        assert(ref_len==alt_len)
    except AssertionError:
        print ref_len,alt_len,index
    return ref_len

fafile = sys.argv[1]
windsize = int(sys.argv[2])
outfile = sys.argv[3]
modelfile = sys.argv[4]
featuredir = sys.argv[5]
featurecode = sys.argv[6]
order = sys.argv[7].split('_')
cpg_range = int(sys.argv[8])

cnt = 0
ref = []
alt = []
ref_anno = []
alt_anno = []
with open(fafile) as f:
    for x in f:
        cnt = (cnt+1)%4
        if cnt == 2:
            ref.append(x.strip())
        if cnt == 0:
            alt.append(x.strip())
        if cnt == 1:
            ref_anno.append(x.strip())
        if cnt == 3:
            alt_anno.append(x.strip())

data = [[ref[idx],alt[idx]] for idx in range(len(ref))]
anno = [[ref_anno[idx],alt_anno[idx]] for idx in range(len(ref))]
print 'Finish preprocess'

cwd = dirname(realpath(__file__))
variant_num = len(data)
convert_dir = join(dirname(fafile),'variant_info')
if not exists(convert_dir):
    makedirs(convert_dir)

if '1' in order:
    outdir = mkdtemp()

    args = [[data[idx],windsize,idx,outdir,anno[idx],cpg_range]for idx in range(variant_num)]
    pool = mp.Pool(processes=20)
    alllength = pool.map(slave,args)
    pool.close()
    pool.join()

    print 'Finish getting CpGs in parallel'
    with open(join(convert_dir,'alllength'),'w') as f:
        for x in alllength:
            f.write('%d\n' %x)
    for cls in ['ref','alt']:
        t_seq = join(convert_dir,'all'+cls)
        t_target = join(convert_dir,'all'+cls+'.t')
        t_h5 = join(convert_dir,'all'+cls+'.h5')
        system('touch ' + t_seq)
        print 'Num of files to concatenate:',str(variant_num)
        print 'Under:',outdir
        with open(t_seq,'w') as f:
            for idx in range(variant_num):
                with open(join(outdir,'test.'+cls+'evl.'+str(idx))) as tf:
                    for x in tf:
                        f.write(x)
        print 'Finish concatenating'

        with open(t_target,'w') as f:
            for idx in range(sum(alllength)):
                f.write('1\n')

        cmd = ' '.join(['python',join(cwd,'helper/embedH5.py'),t_seq,t_target,t_h5,'-b 100000'])
        system(cmd)
    system('rm -r '+outdir)

if '2' in order:

    with open(join(convert_dir,'alllength')) as f:
        alllength = [int(x) for x in f]

    for cls in ['ref','alt']:
        t_h5 = join(convert_dir,'all'+cls+'.h5')
        t_outfile = join(convert_dir,featurecode+'.pred.all'+cls+'.h5')
        cmd = ' '.join(['python',join(cwd,'main.py'),' -p -i',t_h5+'.batch','-s',str(2*windsize+1),'-d',featuredir,'-m',modelfile,'-o',t_outfile])
        print cmd
        system(cmd)

    with open(join(convert_dir,featurecode+'.pred.allref.h5')) as f:
        ref = np.asarray([float(x.strip())for x in f])

    with open(join(convert_dir,featurecode+'.pred.allalt.h5')) as f:
        alt = np.asarray([float(x.strip())for x in f])

    offset = 0
    with open(outfile,'w') as f:
        f.write(featurecode+'\n')
        for idx in range(variant_num):
            new_offset = offset + alllength[idx]
            if new_offset == offset:
                t_ref = [0.001]
                t_alt = [0.001]
            else:
                t_ref = ref[offset:new_offset]
                t_alt = alt[offset:new_offset]
            f.write('%s\n' % ' '.join(map(str,[sum(t_ref),np.mean(t_ref),max(t_ref)])))
            f.write('%s\n' % ' '.join(map(str,[sum(t_alt),np.mean(t_alt),max(t_alt)])))
            offset = new_offset
