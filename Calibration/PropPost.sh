qsub -r y -j y <<EOF

#\$ -V
#\$ -q warp.q@iru3
#\$ -cwd
#\$ -o outputPost.dat
#\$ -e errorPost.dat
#\$ -N "AG1KSamples"
#\$ -pe smp 36

echo "Strating job: AnderNasa"

/apps/MATLAB/R2018b/bin/matlab -nodisplay -r PropagatePosterior

echo "Done with job: AnderModelUp"

EOF

