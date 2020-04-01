qsub -r y -j y <<EOF

#\$ -V
#\$ -q warp.q@iru2
#\$ -cwd
#\$ -o output.dat
#\$ -e error.dat
#\$ -N "AG1KSamples"
#\$ -pe smp 24

echo "Strating job: AnderNasa"

/apps/MATLAB/R2018b/bin/matlab -nodisplay -r runSIRUQ

echo "Done with job: AnderModelUp"

EOF

