java -cp ~/ftp-cbcb/pub/data/PBcR/closure_paper/wgs-package/tools/:. simulateSequencing lambda.fasta 30000 0.0222 0.0075 0.0003 lambda.lens 100  > lambda.5error.fasta
/fs/wrenhomes/sergek/validation/MUMmer3.23/nucmer --maxmatch -l 8 lambda.fasta lambda.5error.fasta
/fs/wrenhomes/sergek/validation/MUMmer3.23/show-coords -lrcTH out.delta |awk '{if ($NF != $(NF-1)) print $0}'|sort -nk12 -nk13 > lambda.5error.vsRef
/cbcb/personal-scratch/sergek/smrtanalysis/current/analysis/bin/blasr lambda.5error.fasta lambda.5error.fasta -bestn 50 -minMatch 2 -maxLCPLength 15 -m 4 -minPctIdentity 50 -maxScore 10000 |awk '{if ($7-$6 > 4000 && $1 > $2) print $0}'|sort -nk1 -nk2 > lambda.5error.vsSelf

java -cp ~/ftp-cbcb/pub/data/PBcR/closure_paper/wgs-package/tools/:. simulateSequencing lambda.fasta 30000 0.037 0.0125 0.0005 lambda.lens 100  > lambda.10error.fasta
/fs/wrenhomes/sergek/validation/MUMmer3.23/nucmer --maxmatch -l 8 lambda.fasta lambda.10error.fasta 
/fs/wrenhomes/sergek/validation/MUMmer3.23/show-coords -lrcTH out.delta |awk '{if ($NF != $(NF-1)) print $0}'|sort -nk12 -nk13 > lambda.10error.vsRef
/cbcb/personal-scratch/sergek/smrtanalysis/current/analysis/bin/blasr lambda.10error.fasta lambda.10error.fasta -bestn 50 -minMatch 2 -maxLCPLength 15 -m 4 -minPctIdentity 50 -maxScore 10000 |awk '{if ($7-$6 > 4000 && $1 > $2) print $0}'|sort -nk1 -nk2 > lambda.10error.vsSelf

java -cp ~/ftp-cbcb/pub/data/PBcR/closure_paper/wgs-package/tools/:. simulateSequencing lambda.fasta 30000 0.066 0.0225 0.0009 lambda.lens 100  > lambda.15error.fasta
/fs/wrenhomes/sergek/validation/MUMmer3.23/nucmer --maxmatch -l 8 lambda.fasta lambda.15error.fasta
/fs/wrenhomes/sergek/validation/MUMmer3.23/show-coords -lrcTH out.delta |awk '{if ($NF != $(NF-1)) print $0}'|sort -nk12 -nk13 > lambda.15error.vsRef
/cbcb/personal-scratch/sergek/smrtanalysis/current/analysis/bin/blasr lambda.15error.fasta lambda.15error.fasta -bestn 50 -minMatch 2 -maxLCPLength 15 -m 4 -minPctIdentity 50 -maxScore 10000 |awk '{if ($7-$6 > 4000 && $1 > $2) print $0}'|sort -nk1 -nk2 > lambda.15error.vsSelf

java -cp ~/ftp-cbcb/pub/data/PBcR/closure_paper/wgs-package/tools/:. simulateSequencing lambda.fasta 30000 0.0814 0.0275 0.0011 lambda.lens 100  > lambda.20error.fasta
/fs/wrenhomes/sergek/validation/MUMmer3.23/nucmer --maxmatch -l 8 lambda.fasta lambda.20error.fasta
/fs/wrenhomes/sergek/validation/MUMmer3.23/show-coords -lrcTH out.delta |awk '{if ($NF != $(NF-1)) print $0}'|sort -nk12 -nk13 > lambda.20error.vsRef
/cbcb/personal-scratch/sergek/smrtanalysis/current/analysis/bin/blasr lambda.20error.fasta lambda.20error.fasta -bestn 50 -minMatch 2 -maxLCPLength 15 -m 4 -minPctIdentity 50 -maxScore 10000 |awk '{if ($7-$6 > 4000 && $1 > $2) print $0}'|sort -nk1 -nk2 > lambda.20error.vsSelf

java -cp ~/ftp-cbcb/pub/data/PBcR/closure_paper/wgs-package/tools/:. simulateSequencing lambda.fasta 30000 0.111 0.0375 0.0015 lambda.lens 100  > lambda.25error.fasta
/fs/wrenhomes/sergek/validation/MUMmer3.23/nucmer --maxmatch -l 8 lambda.fasta lambda.25error.fasta
/fs/wrenhomes/sergek/validation/MUMmer3.23/show-coords -lrcTH out.delta |awk '{if ($NF != $(NF-1)) print $0}'|sort -nk12 -nk13 > lambda.25error.vsRef
/cbcb/personal-scratch/sergek/smrtanalysis/current/analysis/bin/blasr lambda.25error.fasta lambda.25error.fasta -bestn 50 -minMatch 2 -maxLCPLength 15 -m 4 -minPctIdentity 50 -maxScore 10000 |awk '{if ($7-$6 > 4000 && $1 > $2) print $0}'|sort -nk1 -nk2 > lambda.25error.vsSelf

java -cp ~/ftp-cbcb/pub/data/PBcR/closure_paper/wgs-package/tools/:. simulateSequencing lambda.fasta 30000 0.1332 0.045 0.0018 lambda.lens 100  > lambda.30error.fasta
/fs/wrenhomes/sergek/validation/MUMmer3.23/nucmer --maxmatch -l 8 lambda.fasta lambda.30error.fasta
/fs/wrenhomes/sergek/validation/MUMmer3.23/show-coords -lrcTH out.delta |awk '{if ($NF != $(NF-1)) print $0}'|sort -nk12 -nk13 > lambda.30error.vsRef
/cbcb/personal-scratch/sergek/smrtanalysis/current/analysis/bin/blasr lambda.30error.fasta lambda.30error.fasta -bestn 50 -minMatch 2 -maxLCPLength 15 -m 4 -minPctIdentity 50 -maxScore 10000 |awk '{if ($7-$6 > 4000 && $1 > $2) print $0}'|sort -nk1 -nk2 > lambda.30error.vsSelf
