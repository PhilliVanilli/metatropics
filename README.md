# metatropics
a python pipeline for nanopore metagenomics v2.0

git clone --recursive https://github.com/PhilliVanilli/metatropics.git

cd metatropics

conda env create -f meta_dev.yml

conda activate meta

cd jvarkit 
./gradlew sam4weblogo

change threads and demultiplex parameters in code depending on gpu/cpu of the machine

copy host_genomes install dorado and download all models through cmd line
curl "https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.9.1-linux-x64.tar.gz" -o dorado-0.9.1-linux-x64.tar.gz
tar -xzf dorado-0.9.1-linux-x64.tar.gz
dorado-0.9.1-linux-x64/bin/dorado --version

dorado-0.9.1-linux-x64/bin/dorado download --model

copy custom fastq and custom arr toml


copying tar file and unpacking it through remote access doesn't work, one needs to unpack on the machine itself