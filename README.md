# CpGenie
A deep learning method for predicting DNA methylation level of CpG sites from the sequence context, and predicting non-coding variants' effects on DNA methylation.

## Dependencies
+ [Docker](https://www.docker.com/)
+ [NVIDIA-docker](https://github.com/NVIDIA/nvidia-docker)
+ [NVIDIA CUDA](https://developer.nvidia.com/cuda-zone): currently we support CUDA 7.0 and CUDA 8.0.

## Predict DNA methylation level of CpG sites
+ Prepare a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file of the 1001 bp sequence context centered at the CpG you wish to predict, one sequence for one CpG. The 501 nucleotide should be the 'C' of the CpG.
+ Process the FASTA file and make predictions with 50 CpGenie models trained on RRBS datasets from ENCODE immortal cell lines.
	
	```bash
		docker pull haoyangz/cpgenie:CUDA_VER
		mkdir -p OUTPUT_DIR
		nvidia-docker run -u $(id -u) -v FASTA_FILE:/in.fa -v OUTPUT_DIR:/outdir \
				--rm haoyangz/cpgenie:CUDA_VER python main.py ORDER -cpg_fa /in.fa -cpg_out /outdir
	```
	+ `CUDA_VER`: 'cuda7.0' or 'cuda8.0' depending on your NVIDIA driver version.
	+ `FASTA_FILE`: the *absolute path* to the FASTA file of sequences to predict.
	+ `OUTPUT_DIR`: the *absolute path* to the output directory, under which the prediction from each of the 50 CpGenie models will be saved. 
	+ `ORDER`: the following orders can be both used and seperated by space:
		+ `-embed`: data preprocessing. The output will be saved under $OUTPUT_DIR$/embedded.h5.
		+ `-cpg`: make prediction and save under **$OUTPUT_DIR$/CpGenie_pred**.
	+ If you wish to only predict for a subset of the 50 CpGenie models
		+ Download the 50 CpGenie models to a customized directory (e.g. YOUR_MODEL_DIR, an **absolute path**)
			
			```bash
				mkdir -p YOUR_MODEL_DIR
				cd YOUR_MODEL_DIR
				wget http://gerv.csail.mit.edu/CpGenie_models.tar.gz
				tar -zxvf CpGenie_models.tar.gz
			```
		+ Keep only the models you want
		+ Run the following instead:
			
			```bash
				docker pull haoyangz/cpgenie:CUDA_VER
				mkdir -p OUTPUT_DIR
				nvidia-docker run -u $(id -u) -v FASTA_FILE:/in.fa -v OUTPUT_DIR:/outdir -v YOUR_MODEL_DIR:/modeldir \
						--rm haoyangz/cpgenie:CUDA_VER python main.py ORDER -cpg_fa /in.fa -cpg_out /outdir \
							-modeltop /modeldir/models
			```

## Predict the functional score of sequence variants
+ Prepare the variants to score in [VCF](http://www.1000genomes.org/wiki/Analysis/vcf4.0/) format.
+ Score each variant by the predicted impact on DNA methylation in 50 RRBS datasets.

	```bash
		docker pull haoyangz/cpgenie:CUDA_VER
		mkdir -p OUTPUT_DIR
	 	nvidia-docker run -u $(id -u) -v VCF_FILE:/in.vcf -v OUTPUT_DIR:/outdir \
				--rm haoyangz/cpgenie:CUDA_VER python main.py ORDER -var_vcf /in.vcf -var_outdir /outdir
	```
	+ `CUDA_VER`: 'cuda7.0' or 'cuda8.0' depending on your NVIDIA driver version.
	+ `VCF_FILE`: the *absolute path* to the VCF file to score.
	+ `OUTPUT_DIR`: the *absolute path* to the output directory.
	+ `ORDER`: the following orders can be both used and seperated by space:
		+ `-var_prep`: find all CpG sites within 500 bp to each variant and prepare the right format under $OUTPUT_DIR$/CpGenie_processed.
		+ `-var_score`: make predictions on the data preprocessed in the previous step. For each variant, the predicted absolute change of the following metrics are generated for each of the 50 RRBS datasets, resulting in a 250-dim feature vector. The output will be saved as **$OUTPUT_DIR$/CpGenie_var_pred**, where the first line is a header of feature names and then one line for each variant's 250-dim feature vector.
			+ sum of methylation within 500 bp
			+ max methylation within 500 bp 
			+ log odds of the max methylation within 500 bp
			+ mean methylation within 500 bp
			+ log odds of the mean methylation within 500 bp

## Running with CPU
You can run CpGenie on CPU  by replacing all the `nvidia-docker` with `docker`, but it will be extremely slow and thus highly not recommended.
