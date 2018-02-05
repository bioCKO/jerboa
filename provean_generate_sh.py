import sys
# argv[1]: path
# argv[2]: number of genes in each job
# argv[3]: total job numbers

for i in range(1, int(sys.argv[3]) + 1):
    print(i)
    f = open(sys.argv[1] + '/' + 'provean.%s.%d.sh' % (sys.argv[2], i), 'w')
    f.write("""#!/bin/bash
#SBATCH --job-name="{i}_provean"
#SBATCH --output="provean.{i}.out"
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --export=ALL
#SBATCH -t 47:00:00

python /home/hanqing/src/provean.py {num} {i}
""".format(i=i, num=sys.argv[2]))
    f.close()

f = open(sys.argv[1] + '/' + 'total_submit.sh', 'w')
f.write("""
for FILE in provean.*.sh; do
echo ${FILE}
sbatch ${FILE}
sleep 2 # pause to be kind to the scheduler
done
""")
f.close()