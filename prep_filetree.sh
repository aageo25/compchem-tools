# Collect what is needed
# the main folder of the project
echo "Where should I create the file tree?"
read projdir
# setting up the environment
echo "Please give a name to this project:"
read env
mainfolder="$projdir/$env"
# setup the database
echo "Please provide the database file name:"
read db

echo "The main project folder will be: $mainfolder"
echo "The main database is: $db"
echo "Proceed?(yes/no)"
read answer

if [ $answer != 'yes' ]; then
    echo "Aborting"
else
if [ -d "$mainfolder" ]; then
    echo "$mainfolder already exists. Aborting"
else
# creating the folders
mkdir $mainfolder $mainfolder/{dbs,scripts,settings,run,backup}
# creating the source file
cat << EOF > $mainfolder/settings/prepare_env.sh
# setting up the environment
env="$env"
# the main folder of the project
mainfolder="$mainfolder"
# setup the database
structuresdb="$mainfolder/dbs/$db"
# adjust PS1 to avoid mistakes
PS1="\n(\\\$env)\n\033[01;32m\\\$PWD\n${debian_chroot:+($debian_chroot)}\[\033[01;36m\]\u@\h\[\033[00m\]>> "

# source the necessary modules
ml purge > /dev/null 2>&1 # Throw away the warnings from purge
ml iccifort/2019.5.281 impi/2019.7.217 # librearies for VASP
ml GCCcore/8.3.0 # library for Python and tkinter
ml Tkinter/3.7.4 # for making ase gui works
ml VASP/5.4.4-18Apr17-p01-hpc2n

# load python virtual environment
source $HOME/Public/vpyenv/bin/activate

# setting up the calculator
export VASP_PP_PATH="$HOME/calcs/env/vasp_pp"
export VASP_COMMAND="srun vasp_std"
EOF
# creating the run.sh file
cat << EOF > $mainfolder/scripts/run.sh
#!/bin/bash -l
# The -l above is required to get the full environment with modules

#SBATCH -A snic2020-1-41
#SBATCH --mail-type=ALL
#SBATCH -t 0-01:00:00

#SBATCH -J id_xx
#SBATCH -o dat.out

# number of nodes
#SBATCH -n 28  
#SBATCH --ntasks-per-node 28

# prepare bash
source $mainfolder/settings/prepare_env.sh

# RUNNING THE PYTHON CODE
python -u ../run.py \$structuresdb
EOF
# creating worker.sh
cat << 'EOF' > $mainfolder/scripts/worker.sh
#!/bin/bash

calc_id="$(find . -type d | sort -n |awk -F  "/" '{print $2}')"
#calc_id="$(cat list_id.txt)"
#calc_id="$(seq 38 48)"
#calc_id="41 42"

home_dir=$(pwd)

for i in $calc_id ; do
    work_id=${i}
    cd $work_id
    echo ========
    echo $work_id
    echo ========
    if [ -f "WAVECAR" ]; then
        rm WAVECAR # necessary if you restart by hand
    fi
    # submit the calculation
    sbatch --job-name=id_$work_id.ni_g14_alloys run.sh > JobID ; more JobID
    # run a script
#    ase db $structuresdb id=$work_id -k converged=True
#    python ../update_db.py $structuresdb
    cd $home_dir
done
EOF
fi
fi
