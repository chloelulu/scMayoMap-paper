#!/bin/bash

## For internal DB and save the run time
cd /research/bsi/projects/staff_analysis/m216453/scmapping/result/GenomeBiology2019_Tamin/PBMCbench/annotation_res/TabulaMuris/SCSA/
module load python/3.6.5
echo "File,RealTime" > runtimes.csv
for file in *.csv;
do
  real_time=$( { time python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o output2/${file}.txt -m txt ; } 2>&1 | grep real | awk '{print $2}' )
  echo "${file},${real_time}" >> runtimes.csv
done



cd /research/bsi/projects/staff_analysis/m216453/scmapping/result/GenomeBiology2019_Tamin/PBMCbench
module load python/3.6.5
echo "File,RealTime" > runtimes.csv
for file in in *markers*;
do
  real_time=$( { time python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o /research/bsi/projects/staff_analysis/m216453/scmapping/result/GenomeBiology2019_Tamin/PBMCbench/annotation_res/TabulaMuris/SCSA/output2/${file}.txt -m txt ; } 2>&1 | grep real | awk '{print $2}' )
  echo "${file},${real_time}" >> runtimes.csv
done


## For use scMayoDB
cd /research/bsi/projects/staff_analysis/m216453/scmapping/result/GenomeBiology2019_Tamin/PBMCbench/annotation_res/TabulaMuris/SCSA/
module load python/3.6.5
for file in BaronHuman_Pancreas.csv BaronMouse_Pancreas.csv Xin_Pancreas.csv Muraro_Pancreas.csv Segerstolpe_Pancreas.csv Pancreas_smartseq2.csv;
do
  python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o scMayoDB/${file}.txt -m txt -M /research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/pancreas.tsv -N -b
done

for file in Brain_Myeloid_smartseq2.csv Brain_Non-Myeloid_smartseq2.csv;
do
  python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o scMayoDB/${file}.txt -m txt -M /research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/brain.tsv -N -b
done

for file in Heart_smartseq2.csv Heart_and_Aorta_droplet.csv;
do
  python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o scMayoDB/${file}.txt -m txt -M /research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/heart.tsv -N -b
done

for file in Kidney_droplet.csv Kidney_smartseq2.csv;
do
  python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o scMayoDB/${file}.txt -m txt -M /research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/kidney.tsv -N -b
done

for file in Liver_smartseq2.csv Liver_droplet.csv;
do
  python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o scMayoDB/${file}.txt -m txt -M /research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/liver.tsv -N -b
done

for file in Limb_Muscle_droplet.csv Limb_Muscle_smartseq2.csv;
do
  python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o scMayoDB/${file}.txt -m txt -M /research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/muscle.tsv -N -b
done

for file in Lung_droplet.csv Lung_smartseq2.csv;
do
  python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o scMayoDB/${file}.txt -m txt -M /research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/lung.tsv -N -b
done

for file in Large_Intestine_smartseq2.csv;
do
  python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o scMayoDB/${file}.txt -m txt -M /research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/gastrointestinaltract.tsv -N -b
done

for file in Skin_smartseq2.csv;
do
  python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o scMayoDB/${file}.txt -m txt -M /research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/skin.tsv -N -b
done

for file in Marrow_droplet.csv Marrow_smartseq2.csv;
do
  python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o scMayoDB/${file}.txt -m txt -M /research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/bonemarrow.tsv -N -b
done


for file in Spleen_smartseq2.csv Spleen_droplet.csv;
do
  python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o scMayoDB/${file}.txt -m txt -M /research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/spleen.tsv -N -b
done

for file in Thymus_droplet.csv Thymus_smartseq2.csv;
do
  python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o scMayoDB/${file}.txt -m txt -M /research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/thymus.tsv -N -b
done

for file in Bladder_droplet.csv Bladder_smartseq2.csv;
do
  python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o scMayoDB/${file}.txt -m txt -M /research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/bladder.tsv -N -b
done

for file in Mammary_Gland_smartseq2.csv Mammary_Gland_droplet.csv;
do
  python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o scMayoDB/${file}.txt -m txt -M /research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/mammarygland.tsv -N -b
done

for file in Fat_smartseq2.csv;
do
  python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o scMayoDB/${file}.txt -m txt -M /research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/adiposetissue.tsv -N -b
done

for file in Tongue_droplet.csv Tongue_smartseq2.csv;
do
  python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o scMayoDB/${file}.txt -m txt -M /research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/scMayoMapDatabase.tsv -N -b
done

for file in Trachea_smartseq2.csv Trachea_droplet.csv;
do
  python3 /research/bsi/projects/staff_analysis/m216453/SCSA/SCSA.py -d /research/bsi/projects/staff_analysis/m216453/SCSA/whole.db -i ${file} -s seurat -k All -E -g Human -o scMayoDB/${file}.txt -m txt -M /research/bsi/projects/staff_analysis/m216453/SCSA/scMayoDB/lung.tsv -N -b
done
