echo "copying"
cp /root/MorganFP_prediction/entropy_based_datasets/Spectre_2d1d_canonical_smiles_cleaned_v2.zip  /workspace/

echo "start to unzip..."

# unzip -q /root/MorganFP_prediction/entropy_based_datasets/combined_two_datasetes_and_rankingsets_without_empty_oned_NMR.zip -d /workspace/
# unzip -q /workspace/combined_two_datasetes_and_rankingsets_without_empty_oned_NMR.zip -d /workspace/
# unzip -q /workspace/weird_H_and_tautomer_cleaned.zip -d /workspace/
# unzip -q /workspace/RemoveIsomericInfoSMILES.zip -d /workspace/
unzip -q /workspace/Spectre_2d1d_canonical_smiles_cleaned_v2.zip -d /workspace/
# 
# SMART_2D_combined_by_canonical_smiles.zip
# OneD_Only_Dataset.zip


# unzip -q /root/MorganFP_prediction/cleaned_dataset/ranking_sets_with_all_info_molecules.zip -d /workspace

echo "unzip done"


# zip -r combined.zip folder1 folder2
