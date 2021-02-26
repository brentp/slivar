set -e
<<DONE
python plot-final-genome.py slivar-summ-results/genome.sarcoma.gatk.slivar.summary.txt slivar-summ-results/genome.sarcoma.gatk.ch.slivar.summary.txt slivar-summ-results/genome.sarcoma.dv.slivar.summary.txt slivar-summ-results/genome.sarcoma.dv.ch.slivar.summary.txt  slivar-summ-results/genome.ped 2>err
mv figure5-genome-counts.png supp-figure-sarcoma-genome-counts.png
mv figure5-genome-counts.eps supp-figure-sarcoma-genome-counts.eps
#

DONE

python plot-final-genome.py slivar-summ-results/genome.sarcoma.gatk.slivar.impactful.summary.txt slivar-summ-results/genome.sarcoma.gatk.ch.slivar.impactful.summary.txt slivar-summ-results/genome.sarcoma.dv.slivar.impactful.summary.txt slivar-summ-results/genome.sarcoma.dv.ch.slivar.impactful.summary.txt  slivar-summ-results/genome.ped 2>err
mv figure5-genome-counts.png supp-impactful-figure-sarcoma-genome-counts.png
mv figure5-genome-counts.eps supp-impactful-figure-sarcoma-genome-counts.eps

exit

<<DONE

python plot-final-genome.py slivar-summ-results/genome.sarcoma.gatk.slivar.genic.summary.txt slivar-summ-results/genome.sarcoma.gatk.ch.slivar.genic.summary.txt slivar-summ-results/genome.sarcoma.dv.slivar.genic.summary.txt slivar-summ-results/genome.sarcoma.dv.ch.slivar.genic.summary.txt  slivar-summ-results/genome.ped 2>err
mv figure5-genome-counts.png supp-genic-figure-sarcoma-genome-counts.png
mv figure5-genome-counts.eps supp-genic-figure-sarcoma-genome-counts.eps



python plot-final-genome.py slivar-summ-results/genome.RGP.gatk.slivar.summary.txt slivar-summ-results/genome.RGP.gatk.ch.slivar.summary.txt slivar-summ-results/genome.RGP.dv.slivar.summary.txt slivar-summ-results/genome.RGP.dv.ch.slivar.summary.txt  slivar-summ-results/genome.ped 2>err
mv figure5-genome-counts.png supp-figure-RGP-genome-counts.png
mv figure5-genome-counts.eps supp-figure-RGP-genome-counts.eps

DONE

python plot-final-genome.py slivar-summ-results/genome.RGP.gatk.slivar.impactful.summary.txt slivar-summ-results/genome.RGP.gatk.ch.slivar.impactful.summary.txt slivar-summ-results/genome.RGP.dv.slivar.impactful.summary.txt slivar-summ-results/genome.RGP.dv.ch.slivar.impactful.summary.txt  slivar-summ-results/genome.ped 2>err
mv figure5-genome-counts.png supp-impactful-figure-RGP-genome-counts.png
mv figure5-genome-counts.eps supp-impactful-figure-RGP-genome-counts.eps


<<DONE

python plot-final-genome.py slivar-summ-results/genome.RGP.gatk.slivar.genic.summary.txt slivar-summ-results/genome.RGP.gatk.ch.slivar.genic.summary.txt slivar-summ-results/genome.RGP.dv.slivar.genic.summary.txt slivar-summ-results/genome.RGP.dv.ch.slivar.genic.summary.txt  slivar-summ-results/genome.ped 2>err
mv figure5-genome-counts.png supp-genic-figure-RGP-genome-counts.png
mv figure5-genome-counts.eps supp-genic-figure-RGP-genome-counts.eps
#
DONE
