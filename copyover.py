from os import system,makedirs
from os.path import join,exists

model = 'seq_128x3_5_5_2f_simple'
with open('log') as f:
    for x in f:
        expt = x.split()[-1]+'.combined.bed.1001bp'
        outdir = join('models',expt,model)
        if not exists(outdir):
            makedirs(outdir)
        system(' '.join(['cp',join('../methy/FEATURE1/',expt,model,'*_best_archit.json'),outdir]))
        system(' '.join(['cp',join('../methy/FEATURE1/',expt,model,'*_bestmodel_weights.h5'),outdir]))
