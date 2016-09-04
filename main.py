from __future__ import print_function
import theano
import numpy as np,sys,h5py,cPickle,argparse,subprocess
from hyperopt import Trials, STATUS_OK, tpe
from hyperas import optim
from os.path import join,dirname,basename,exists,realpath,isdir
from os import system,chdir,getcwd,listdir,makedirs
from keras.models import model_from_json
from keras.optimizers import SGD
from tempfile import mkdtemp,NamedTemporaryFile

cwd = dirname(realpath(__file__))
t_wind = 500
genomefile = join(cwd,'data/hg19.in')
sizefile = join(cwd,'data/hg19.size')
variant_script = join(cwd,'variant.py')

def parse_args():
    parser = argparse.ArgumentParser(description="Launch a list of commands on EC2.")
    parser.add_argument("-y", "--hyper", dest="hyper", default=False, action='store_true',help="")
    parser.add_argument("-t", "--train", dest="train", default=False, action='store_true',help="")
    parser.add_argument("-e", "--eval", dest="eval", default=False, action='store_true',help="")
    parser.add_argument("-p", "--predit", dest="predict", default=False, action='store_true',help="")
    parser.add_argument("-i", "--infile", dest="infile", default='',help="")
    parser.add_argument("-d", "--topdir", dest="topdir",default=cwd,help="")
    parser.add_argument("-s", "--datasize", dest="datasize",default=1001,help="")
    parser.add_argument("-c", "--datacode", dest="datacode",default='data',help="")
    parser.add_argument("-v", "--cvfold", dest="cvfold",default='0',help="")
    parser.add_argument("-m", "--model", dest="model",default=join(cwd,'cnn/seq_128x3_5_5_2f_simple'),help="")
    parser.add_argument("-o", "--outfile", dest="outfile",default='',help="")
    parser.add_argument("-x", "--prefix", dest="prefix",default='',help="")
    parser.add_argument("-hi", "--hyperiter", dest="hyperiter",default=9,help="")
    parser.add_argument("-cpg", dest="cpg",default=False,action='store_true',help="")
    parser.add_argument("-embed", dest="embed",default=False,action='store_true',help="")
    parser.add_argument("-var_prep", dest="var_prep",default=False,action='store_true',help="")
    parser.add_argument("-var_score", dest="var_score",default=False,action='store_true',help="")
    parser.add_argument("-modeltop", dest="modeltop",default=join(cwd,'models'),help="")
    parser.add_argument("-cpg_fa", "--cpg_fa",dest="cpg_fa",default='',help="")
    parser.add_argument("-cpg_out", "--cpg_out",dest="cpg_out",default='',help="")
    parser.add_argument("-var_vcf", "--var_vcf",dest="var_vcf",default='',help="")
    parser.add_argument("-var_vcf_tmp", "--var_vcf_tmp",dest="var_vcf_tmp",default='',help="")
    parser.add_argument("-var_outdir", "--var_outdir",dest="var_outdir",default='',help="")
    return parser.parse_args()

def createdir(mydir):
    if exists(mydir):
        system('rm -r ' + mydir)
    makedirs(mydir)

if __name__ == "__main__":
    args = parse_args()
    #system(' '.join(['cp',args.model+'.template','.']))
    #topdir = join(args.topdir,'CV'+args.cvfold)
    topdir = args.topdir
    DATASIZE = args.datasize
    model_arch = basename(args.model)
    data_code = args.datacode

    outdir = join(topdir,model_arch)
    if args.hyper or args.train:
        if not exists(outdir):
            makedirs(outdir)
    architecture_file = join(outdir,model_arch+'_best_archit.json')
    optimizer_file = join(outdir,model_arch+'_best_optimer.pkl')
    weight_file = join(outdir,model_arch+'_bestmodel_weights.h5')
    data1prefix = join(topdir,data_code+args.prefix)
    evalout = join(outdir,model_arch+'_eval.txt')

    tmpdir = mkdtemp()
    with open(args.model+'.template') as f,open(join(tmpdir,model_arch+'.py'),'w') as fout:
        for x in f:
            newline = x.replace('DATACODE',data_code)
            newline = newline.replace('TOPDIR',topdir)
            newline = newline.replace('DATASIZE',str(DATASIZE))
            newline = newline.replace('MODEL_ARCH',model_arch)
            newline = newline.replace('PREFIX',args.prefix)
            fout.write(newline)

    sys.path.append(tmpdir)
    mymodel = __import__(model_arch)
    allmodels = [x for x in listdir(args.modeltop) if isdir(join(args.modeltop,x))]

    if args.hyper:
        ## Hyper-parameter tuning
        MAX_EVAL = args.hyperiter
        best_run, best_model = optim.minimize(model=mymodel.model,data=mymodel.data,algo=tpe.suggest,max_evals=MAX_EVAL,trials=Trials())
        best_archit,best_optim = best_model
        open(architecture_file, 'w').write(best_archit)
        cPickle.dump(best_optim,open(optimizer_file,'wb') )

    if args.train:
        ### Training
        from keras.models import model_from_json
        from keras.callbacks import ModelCheckpoint

        model = model_from_json(open(architecture_file).read())
        best_optim = cPickle.load(open(optimizer_file,'rb'))
        model.compile(loss='binary_crossentropy', optimizer=best_optim,metrics=['accuracy'])

        checkpointer = ModelCheckpoint(filepath=weight_file, verbose=1, save_best_only=True)
        trainsample = join(topdir,data_code)+'.target*train'
        train_size = int(subprocess.check_output('wc '+trainsample, shell=True).split()[0])
        validsample = join(topdir,data_code)+'.target*valid'
        valid_size = int(subprocess.check_output('wc '+validsample, shell=True).split()[0])
        trainbatch_num = int(subprocess.check_output('ls '+data1prefix+'.train.h5.batch* | wc -l', shell=True).split()[0])
        validbatch_num = int(subprocess.check_output('ls '+data1prefix+'.valid.h5.batch* | wc -l', shell=True).split()[0])
        TRAIN_EPOCH = 20
        BATCHSIZE = 100
        history_callback = model.fit_generator(mymodel.BatchGenerator2(BATCHSIZE,trainbatch_num,'train',topdir,data_code)\
        		    ,train_size,TRAIN_EPOCH,validation_data=mymodel.BatchGenerator2(BATCHSIZE,validbatch_num,'valid',topdir,data_code)\
        			    ,nb_val_samples=valid_size,callbacks = [checkpointer])
        system('touch '+join(outdir,model_arch+'.traindone'))
        loss_history = np.asarray(history_callback.history["loss"])
        acc_history = np.asarray(history_callback.history["acc"])
        val_loss_history = np.asarray(history_callback.history["val_loss"])
        val_acc_history = np.asarray(history_callback.history["val_acc"])
        all_hist = np.hstack((loss_history.reshape(len(loss_history),1),acc_history.reshape(len(acc_history),1),val_loss_history.reshape(len(val_loss_history),1),val_acc_history.reshape(len(val_acc_history),1)))
        np.savetxt(join(outdir,model_arch+".training_history.txt"), all_hist,delimiter = "\t",header='loss\tacc\tval_loss\tval_acc')

    if args.eval:
        ## Evaluate
        model = model_from_json(open(architecture_file).read())
        model.load_weights(weight_file)
        best_optim = cPickle.load(open(optimizer_file,'rb'))
        model.compile(loss='binary_crossentropy', optimizer=best_optim,metrics=['accuracy'])

        pred = np.asarray([])
        y_true = np.asarray([])
        testbatch_num = int(subprocess.check_output('ls '+data1prefix+'.test.h5.batch* | wc -l', shell=True).split()[0])
        for X1_train,Y_train in mymodel.BatchGenerator(testbatch_num,'test',topdir,data_code):
            pred = np.append(pred,[x[0] for x in model.predict(X1_train)])
            y_true = np.append(y_true,[x[0] for x in Y_train])

        from sklearn.metrics import accuracy_score,auc,roc_curve
        fpr, tpr, thresholds = roc_curve(y_true, pred)
        t_auc = auc(fpr, tpr)
        print('Test AUC:',t_auc)

        y_pred = np.zeros(len(pred))
        y_pred[pred>0.5] = 1
        t_acc = accuracy_score(y_true,y_pred)
        print('Test accuracy:',t_acc)

        with open(evalout,'w') as f:
            f.write('%f\n' % t_auc)
            f.write('%f\n' % t_acc)

    if args.predict:
        model = model_from_json(open(architecture_file).read())
        model.load_weights(weight_file)
        model.compile(loss='binary_crossentropy', optimizer=SGD(lr=0.01, momentum=0.9, nesterov=True),metrics=['accuracy'])

        predict_batch_num = int(subprocess.check_output('ls '+args.infile+'* | wc -l', shell=True).split()[0])
        print('Total number of batch to predict:',predict_batch_num)
        if args.outfile == '':
            outfile = join(dirname(args.infile),'pred.'+basename(args.infile))
        else:
            outfile = args.outfile
        if not exists(dirname(outfile)):
            makedirs(dirname(outfile))
        with open(outfile,'w') as f:
            for i in range(predict_batch_num):
                data1f = h5py.File(args.infile+str(i+1),'r')
                data1 = data1f['data']
                pred = model.predict(data1)
                for x in pred:
                    f.write('%f\n' % x[0])
    if args.embed:
        assert(args.cpg_out != '')
        t_outdir = join(args.cpg_out,'embedded.h5')
        if not exists(t_outdir):
            makedirs(t_outdir)
        t_outfile = join(t_outdir,'data')
        tsv_file = NamedTemporaryFile(delete=False).name
        labelfile = NamedTemporaryFile(delete=False).name

        system(' '.join(['paste - - -d\' \'','<',args.cpg_fa,'>',tsv_file]))
        sample_size = int(subprocess.check_output('wc '+tsv_file, shell=True).split()[0])

        with open(labelfile,'w') as f:
            for x in range(sample_size):
                f.write('1\n')

        system(' '.join(['python helper/embedH5.py',tsv_file,labelfile,t_outfile,'-b 100000']))
        system('rm -r '+tsv_file)
        system('rm -r '+labelfile)

    if args.cpg:
        assert(args.cpg_out != '')
        predict_batch_num = int(subprocess.check_output('ls '+ join(args.cpg_out,'embedded.h5','data.batch*')+' | wc -l', shell=True).split()[0])
        print('Total number of batch to predict:',predict_batch_num)
        t_outdir = join(args.cpg_out,'CpGenie_pred')
        createdir(t_outdir)

        for mymodel in allmodels:
            print(mymodel)
            t_modeldir = join(args.modeltop,mymodel,model_arch)
            t_architecture_file = join(t_modeldir,model_arch+'_best_archit.json')
            t_weight_file = join(t_modeldir,model_arch+'_bestmodel_weights.h5')

            model = model_from_json(open(t_architecture_file).read())
            model.load_weights(t_weight_file)
            model.compile(loss='binary_crossentropy', optimizer=SGD(lr=0.01, momentum=0.9, nesterov=True),metrics=['accuracy'])

            with open(join(t_outdir,mymodel),'w') as f:
                for i in range(predict_batch_num):
                    data1f = h5py.File(join(args.cpg_out,'embedded.h5','data.batch'+str(i+1)),'r')
                    data1 = data1f['data']
                    pred = model.predict(data1)
                    for x in pred:
                        f.write('%f\n' % x[0])

    if args.var_prep:
        assert(args.var_vcf != '')
        assert(args.var_outdir != '')
        if args.var_vcf_tmp == '':
            args.var_vcf_tmp = join(args.var_outdir,'CpGenie_processed')
        vcf_chr_dir = join(args.var_vcf_tmp,'vcf_chrsplit')
        fa_split_dir = join(args.var_vcf_tmp,'fasta')

    	### for each file, split again by chr
        createdir(vcf_chr_dir)
        system(' '.join(['python',join(cwd,'helper','splitVCF.py'),args.var_vcf,vcf_chr_dir]))

        ### VCF to FA
        vcf2fasta_script = join(cwd,'helper/vcf2fasta2.R')
        system(' '.join(['Rscript',vcf2fasta_script,vcf_chr_dir,genomefile,sizefile,fa_split_dir,str(t_wind*2),str(t_wind*2),'T']))

        ### Find CpGs around variants
        fafile = join(fa_split_dir,'all.fa')
        system(' '.join(['python',variant_script,fafile,str(t_wind),'1','1','1','1','1']))

    if args.var_score:
        assert(args.var_vcf != '')
        assert(args.var_outdir != '')
        if args.var_vcf_tmp == '':
            args.var_vcf_tmp = join(args.var_outdir,'CpGenie_processed')
        fa_split_dir = join(args.var_vcf_tmp,'fasta')
        fafile = join(fa_split_dir,'all.fa')

        scores_topdir = join(args.var_outdir,'per_allele')
        if not exists(scores_topdir):
            makedirs(scores_topdir)

        for mymodel in allmodels:
            outfile = join(scores_topdir,mymodel)
            featurecode= mymodel + '_' + basename(args.model)
            system(' '.join(['python',variant_script,fafile,str(t_wind),outfile,args.model,join(args.modeltop,mymodel),featurecode,'2']))


        eqtl_score_processed_dir_allencode  = join(args.var_outdir,'CpGenie_var_pred')

        def logfold(ref,alt,pseudo):
            if ref == 1.0 or (1-ref) == 1.0:
                ref += np.sign(0.5 - ref)*pseudo
            if alt == 1.0 or (1-alt) == 1.0:
                alt += np.sign(0.5 - alt)* pseudo
            return abs(np.log(ref/(1-ref))-np.log(alt/(1-alt)))

        pseudo = 1e-6
        def gwas_feature_new(ref,alt):
            return [abs(ref[0]-alt[0]),abs(ref[1]-alt[1]),logfold(ref[1],alt[1],pseudo), abs(ref[2]-alt[2]),logfold(ref[2],alt[2],pseudo)]

        with open(join(eqtl_score_processed_dir_allencode),'w') as fout:
            t_data = None
            for mymodel in allmodels:
                with open(join(scores_topdir,mymodel)) as fin:
                    fin.readline()
                    t_t_data = np.asarray([map(float,x.split()) for x in fin])
                    t_feature = []
                    for x in range(len(t_t_data)/2):
                        ref = t_t_data[2*x]
                        alt = t_t_data[2*x+1]
                        t_feature.append(gwas_feature_new(ref,alt))
                    t_feature = np.asarray(t_feature)
                    t_data = t_feature if t_data is None else np.hstack((t_data,t_feature))
            header = [[mymodel+'_absdiff_sum',mymodel+'_absdiff_mean',mymodel+'_absdiff_LO_mean',mymodel+'_absdiff_max',mymodel+'_absdiff_LO_max'] for mymodel in allmodels]
            fout.write('%s\n' % '\t'.join([y for x in header for y in x]))
            for x in t_data:
                fout.write('%s\n' % '\t'.join(map(str,x)))

    system('rm -r ' + tmpdir)
