# CpGenie
A deep learning method for predicting DNA methylation level of CpG sites from the sequence context, and predicting non-coding variants' effects on DNA methylation.

## Dependencies
+ [Docker](https://www.docker.com/)
+ NVIDIA 346.46 driver

## Predict DNA methylation level of CpG sites
+ Prepare a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file of the 1001 bp sequence context centered at the CpG you wish to predict, one sequence for one CpG. The 501 nucleotide should be the 'C' of the CpG.
+ Embed the FASTA file to HDF5 format readable by MethylDecode and make predictions with MethylDecoder models trained on 50 ENCODE RRBS datasets
	
	```
		docker pull haoyangz/cpgenie
		docker run --device /dev/nvidiactl --device /dev/nvidia-uvm MOREDEVICE \
					-v FASTA_FILE:/in.fa -v OUTPUT_DIR:/outdir haoyangz/cpgenie \
					python main.py ORDER -cpg_fa /in.fa -cpg_out /outdir
	```
	+ `FASTA_FILE`: the *absolute path* to the FASTA file of sequences to predict
	+ `OUTPUT_DIR`: the *absolute path* to the output directory, under which the prediction from each of the 50 MethylDecoder models will be saved. 
	+ `ORDER`: the following orders can be both used and seperated by space
		+ `-embed`: data preprocessing. The output will be saved under $OUTPUT_DIR$/embedded.h5
		+ `-cpg`: make prediction. 
	+ `MOREDEVICE`: For each of the GPU device available on your machine, append one "--device /dev/nvidiaNUM" where NUM is the device index. For hsf1/hsf2 in  Gifford Lab, since there are three GPUs, it should be :
                                                                                                                                                                                                                         
    	```
    	--device /dev/nvidia0 --device /dev/nvidia1 --device /dev/nvidia2
    	```

## Predict the functional score of sequence variants
+ Prepare the variants to score in [VCF](http://www.1000genomes.org/wiki/Analysis/vcf4.0/) format.
+ Score each variant by the predicted impact on DNA methylation in 50 RRBS datasets.

	```
		docker pull haoyangz/cpgenie
	 	docker run --device /dev/nvidiactl --device /dev/nvidia-uvm MOREDEVICE \
	                 -v VCF_FILE:/in.vcf -v OUTPUT_DIR:/outdir haoyangz/cpgenie \
					python main.py ORDER -var_vcf /in.vcf -var_outdir /outdir
	```
	+ `VCF_FILE`: the *absolute path* to the VCF file to score
	+ `OUTPUT_DIR`: the *absolute path* to the output directory
	+ `ORDER`: the following orders can be both used and seperated by space
		+ `-var_prep`: find all CpG sites within 500 bp to each variant and prepare the right format under $OUTPUT_DIR$/CpGenie_processed
		+ `-var_score`: make predictions on the data preprocessed in the previous step. For each variant, the predicted absolute change of the following metrics are generated for each of the 50 RRBS datasets, resulting in a 250-dim feature vector. The output will be saved as $OUTPUT_DIR$/CpGenie_pred, where the first line is a header of feature names and then one line for each variant's 250-dim feature vector.
			+ sum of methylation within 500 bp
			+ max methylation within 500 bp 
			+ log odds of the max methylation within 500 bp
			+ mean methylation within 500 bp
			+ log odds of the mean methylation within 500 bp
	+ `MOREDEVICE`: same as above
