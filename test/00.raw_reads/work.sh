/public/agis/fanwei_group/fanwei/agrocode/novometa/pirs_simple/pirs diploid  -i Ecoli_MG1655_uid57779_genome.fa -c 0  -o Ecoli_MG1655_uid57779_genome.fa


 /public/agis/fanwei_group/fanwei/agrocode/novometa/pirs_simple/pirs simulate -s /public/agis/fanwei_group/fanwei/agrocode/novometa/pirs_simple/Profiles/humNew.PE250.matrix -l 250 -x 20 -m 400  -i Ecoli_MG1655_uid57779_genome.fa -I Ecoli_MG1655_uid57779_genome.fa.snp.indel.invertion.fa -o Ecoli_readlen250_insert400_20X 2> log1 &
  /public/agis/fanwei_group/fanwei/agrocode/novometa/pirs_simple/pirs simulate -s /public/agis/fanwei_group/fanwei/agrocode/novometa/pirs_simple/Profiles/humNew.PE250.matrix -l 250 -x 20 -m 800  -i Ecoli_MG1655_uid57779_genome.fa -I Ecoli_MG1655_uid57779_genome.fa.snp.indel.invertion.fa -o Ecoli_readlen250_insert800_20X  2> log2 &


