# metatropics
a python pipeline for nanopore metagenomics v2.0

git clone --recursive https://github.com/PhilliVanilli/metatropics.git

cd metatropics

conda env create -f meta_dev.yml

conda activate meta

cd jvarkit 
./gradlew sam4weblogo

change threads and demultiple parameters in code depending on gpu/cpu of the machine

copy host_genomes and guppy 


