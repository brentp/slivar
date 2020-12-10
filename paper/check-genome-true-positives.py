import toolshed
import collections
import hashlib

slivar_found = collections.defaultdict(list)

# first file is tsv from slivar

for d in toolshed.reader(1):
    slivar_found[d["sample_id"]].append(d["chr:pos:ref:alt"])

shared_affected = 0
shared_affected_solved = 0
exact_match = 0
sv_del = 0
oocyte = 0
indel_plus_sv_comphet = 0

# 2nd file is participant_details.tsv 
for d in toolshed.reader(2):
    sample_id = d["entity:participant_details_id"]
    if sample_id not in slivar_found: continue
    if d["07_affected_status"] != "Affected": continue
    shared_affected += 1

    key = "chr%s:%s:%s:%s" % (d["13_Chrom-1"], d["14_Pos-1"], d["15_Ref-1"], d["16_Alt-1"])
    if key == "chr:::": continue
    shared_affected_solved += 1

    if key in slivar_found[sample_id]:
        print("OK", sample_id, key)
        exact_match += 1
    else:
        sha = hashlib.sha256(sample_id.encode()).hexdigest()
        #print("WTF", sample_id, key)
        if key.endswith("del"):
            sv_del += 1
        elif sha in (
                "c1b669b32e2b899a15bcd8e5d3e2cc9f5eb457a1b8a1c27fce2ab8f26750c050",
                "8145446cdae4964156aefbb0eb1ab03f2866905de14942cffedc48b782de5086"):
            oocyte += 1
        elif sha in (
               "2b2f722dcb68c22c654c4cc2a9d8db8bda08a88461d5d5d7d89c985ba726eb62",
               "c52f9645ec80ce4c0c712bb4d8fad246e864b04d53ba94e9a29e3aac15f1985c",
                ):
            indel_plus_sv_comphet += 1
        elif sha == "6503b96da638ddab193fa7cbc1e5c20d626f5d2dda7dabe74b09ed4d3e2a677f":
                print("mom:0/0 dad:0/1 kid:1/1 (because of sv deletion but other filters passed)")
        elif sha == "8d853c417e5d02c5362e0ddce3666689f2feb8672e2293ff9483d4bd6d7ead42":
            print(sample_id, key)
            print("X. ref: CACCCTCCACGAT")
            print("X. reported by RGP: pos:802 var:TCCAC/A")
            print("X. found by our pipelines:")
            print("X. pos:802  CCCT/C")
            print("X. pos:808  AC/A")
        else:
            print("BAD", sample_id, key)
            1/0


print("shared_affected", shared_affected)
print("shared_affected_solved: ", shared_affected_solved)
print("exact_match:", exact_match)
print("SV deletion (not sought here):", sv_del)
print("autosome het 2 girls shared with dad", oocyte)
print("comphet missed because 1 side was deletion:", indel_plus_sv_comphet)


"""
entity:participant_details_id	01_project_id	02_family_id	03_Individual_ID	06_sex	07_affected_status	08_phenotype_description	09_hpo_present	29_Date-Uploaded	04_paternal_id	05_maternal_id	11_Gene-1	12_Zygosity-1	13_Chrom-1	14_Pos-1	15_Ref-1	16_Alt-1	17_hgvsc-1	18_hgvsp-1	19_Transcript-1	10_hpo_absent	20_Gene-2	21_Zygosity-2	22_Chrom-2	23_Pos-2	24_Ref-2	25_Alt-2	26_hgvsc-2	27_hgvsp-2	28_Transcript-2	31_Notes	30_Other_seq_data
RGP_1003_3	Rare Genomes Project_Genomes	RGP_1003	RGP_1003_3	Male	Affected	Limb-girdle muscular dystrophy	HP:0003236 (Elevated serum creatine kinase)|HP:0012378 (Fatigue)|HP:0003325 (Limb-girdle muscle weakness)|HP:0003701 (Proximal muscle weakness)	3/1/2019																							
RGP_1004_3	Rare Genomes Project_Genomes	RGP_1004	RGP_1004_3	Female	Affected	Limb-girdle muscular dystrophy	HP:0012432 (Chronic fatigue)|HP:0006785 (Limb-girdle muscular dystrophy)|HP:0001324 (Muscle weakness)|HP:0003202 (Skeletal muscle atrophy)	3/1/2019																							
RGP_1004_4	Rare Genomes Project_Genomes	RGP_1004	RGP_1004_4	Female	Affected	Limb-girdle muscular dystrophy	HP:0006785 (Limb-girdle muscular dystrophy)	3/1/2019																							
RGP_1004_5	Rare Genomes Project_Genomes	RGP_1004	RGP_1004_5	Female	Affected	Limb-girdle muscular dystrophy	HP:0006785 (Limb-girdle muscular dystrophy)	3/1/2019																							
RGP_1006_3	Rare Genomes Project_Genomes	RGP_1006	RGP_1006_3	Male	Affected	Myopathy	HP:0002355 (Difficulty walking)|HP:0003473 (Fatigable weakness)|HP:0002359 (Frequent falls)|HP:0030237 (Hand muscle weakness)|HP:0007340 (Lower limb muscle weakness)|HP:0001324 (Muscle weakness)|HP:0003484 (Upper limb muscle weakness)	3/1/2019																							
RGP_1012_1	Rare Genomes Project_Genomes	RGP_1012	RGP_1012_1	Female	Unaffected	Overgrowth; autism		3/1/2019																							
RGP_1012_2	Rare Genomes Project_Genomes	RGP_1012	RGP_1012_2	Male	Unaffected	Overgrowth; autism		8/16/2019																							
RGP_1012_3	Rare Genomes Project_Genomes	RGP_1012	RGP_1012_3	Male	Affected	Overgrowth; autism	HP:0000729 (Autistic behavior)|HP:0001548 (Overgrowth)	3/1/2019	RGP_1012_2	RGP_1012_1																					
RGP_1013_3	Rare Genomes Project_Genomes	RGP_1013	RGP_1013_3	Male	Affected	Myopathy	HP:0003198 (Myopathy)	3/1/2019																							
"""
