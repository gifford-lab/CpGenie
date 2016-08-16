# MethylDecoder

## Dependencies
+ [Docker](https://www.docker.com/)
+ NVIDIA 346.46 driver

## Predict DNA methylation of single CpG
+ Prepare a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file of the 1001 bp sequence context centered at the CpG you wish to predict, one sequence for one CpG. The 501 nucleotide should be the 'C' of the CpG.
+ Embed the FASTA file to HDF5 format readable by MethylDecode and make predictions with MethylDecoder models trained on 50 ENCODE RRBS datasets
	```
		docker run --device /dev/nvidiactl --device /dev/nvidia-uvm MOREDEVICE -v FASTA_TOPDIR:/indir -v OUTPUT_DIR:/outdir haoyangz/MethylDecoder python main.py ORDER -cpg_fa FASTA_NAME -cpg_out /outdir
	```
	+ `MOREDEVICE`: For each of the GPU device available on your machine,append one "--device /dev/nvidiaNUM" where NUM is the number of that device. For hsf1/hsf2 in  Gifford Lab, since there are three GPU, it should be :

		```
		--device /dev/nvidia0 --device /dev/nvidia1 --device /dev/nvidia2
		```
	+ `FASTA_TOPDIR`: the *absolute path* to the top directory of the FASTA file
	+ `FASTA_NAME`: the filename of the FASTA file
	+ `OUTPUT_DIR`: the *absolute path* to the output directory, under which the prediction from each of the 50 MethylDecoder models will be saved. 
	+ `ORDER`: '-embed' for feature preparation and '-cpg' for prediction. They can be both used and seperated by space. When '-embed' is included, the converted features will be dumped under $FASTA_TOPDIR$/$FASTA_NAME$.embed_h5 for every 100k sequences.

