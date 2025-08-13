conda create -n ipyrad_env ##create ipyrad environment
conda activate ipyrad_env  ##activate env
conda install ipyrad -c conda-forge -c bioconda  ##install ipyrad
##Now get your params file ready "ipyrad_step1.txt"  ###
##Begin Ipyrad from step 1 through till step 7###
ipyrad -p params_step1.txt -s 1
ipyrad -p params_step1.txt -s 2
ipyrad -p params_step1.txt -s 3
ipyrad -p params_step1.txt -s 4
ipyrad -p params_step1.txt -s 5
ipyrad -p params_step1.txt -s 6
ipyrad -p params_step1.txt -s 7
