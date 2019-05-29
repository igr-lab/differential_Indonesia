import pandas as pd

# Open Illumina EPIC probe annotations
probes = pd.read_csv("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/MethylationEPIC_v-1-0_B4.csv", skiprows=7, sep=",")

villages = ["ANKvsWNG",
"ANKvsRIN",
"ANKvsHPM",
"ANKvsPDT",
"WNGvsRIN",
"WNGvsHPM",
"WNGvsPDT",
"RINvsHPM",
"RINvsPDT",
"HPMvsPDT",
"ANKvsMDB",
"ANKvsTLL",
"WNGvsMDB",
"WNGvsTLL",
"RINvsMDB",
"RINvsTLL",
"HPMvsMDB",
"HPMvsTLL",
"PDTvsMDB",
"PDTvsTLL",
"ANKvsMPI",
"WNGvsMPI",
"RINvsMPI",
"HPMvsMPI",
"PDTvsMPI",
"MDBvsMPI",
"TLLvsMPI",
"MDBvsTLL",
"SMBvsMTW",
"SMBvsMPI",
"MTWvsMPI"]
islands = ["SMBvsMTW", "SMBvsMPI", "MTWvsMPI"]
all_comparisons = villages + islands

def annotate_positions(comparison):
    # Open DMP summary statistics as a pandas data frame
    dmp_results = pd.read_csv("/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/DMP_" + comparison + "_adjp1_logFC0_DeconCell_new.txt", sep="\t")
    # Add a new column
    dmp_results["IlmnID"] = dmp_results.index
    dmp_probe_info = probes[probes['IlmnID'].isin(list(dmp_results.index))]
    # Combine by IlmnID
    annotated_dmp_results = pd.merge(dmp_results, dmp_probe_info, on='IlmnID', how='outer')
    # Write to a text file
    path = "/Users/hnatri/Dropbox (ASU)/Indonesian_methylation/DMP_" + comparison + "_adjp1_logFC0_DeconCell_new_annotated_positions.txt"
    pd.DataFrame.to_csv(annotated_dmp_results, path, sep="\t", index=False)

for comparison in all_comparisons:
    annotate_positions(comparison)
