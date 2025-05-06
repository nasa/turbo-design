#/bin/sh

# Setup miniconda 
__conda_setup="$('/home/pjuangph/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/pjuangph/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/pjuangph/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/pjuangph/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

eval "conda activate dev"
eval "/engapps/Tecplot.360/tecplot2024R1/360ex_2024r1/bin/tec360-env -- python post_process.py"

echo "Tecplot post processed complete"