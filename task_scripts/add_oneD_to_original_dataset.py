import os, torch, pickle, json, random
from rdkit import Chem
import tqdm
import tracemalloc
import logging
from collections import defaultdict

# ### START debugging stuff
# # Setup logging
# logging.basicConfig(filename='debug.log', level=logging.DEBUG)

# # Start tracing memory allocations
# tracemalloc.start()

# def check_memory(snapshot, key_type='lineno'):
#     top_stats = snapshot.statistics(key_type)
#     logging.debug("[Top 10]")
#     for stat in top_stats[:10]:
#         logging.debug(stat)
        
# ### END debugging stuff



def get_canonical_smiles(datum):
    """
    Converts a SMILES string to its canonical form.


    Returns:
    str: The canonical SMILES string.
    """
    smiles = datum['smiles']
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise Exception("Invalid SMILES string")
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        return canonical_smiles
    except:
        try:
            inchi = datum['inchi']
            mol = Chem.MolFromInchi(inchi)
            if mol is None:
                raise Exception("Invalid Inchi string")
            canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
            return canonical_smiles
        except:
            # # print("BAD: ", smiles, inchi)
            return None
        

canonical_smiles_to_path = {}
for split in ["test", "train", "val"]:
    with open(f'/workspace/SMILES_dataset/{split}/SMILES/index.pkl', 'rb') as file:
        smiles_mapping = pickle.load(file)
    os.makedirs(f'/workspace/SMILES_dataset/{split}/oneD_NMR', exist_ok=True)
    for num,smiles in smiles_mapping.items():
        canonical_smiles_to_path[smiles] = f'/workspace/SMILES_dataset/{split}/oneD_NMR/{str(num)}.pt'
        
# # print(list(canonical_smiles_to_path.values())[:4])
# exit(0)
        
# construct a dict   NP_to_smiles
# if path exists
if os.path.exists('/workspace/NP_to_smiles_and_names.pkl'):
    # print("load from existing file")
    NP_to_smiles, NP_to_names = pickle.load(open('/workspace/NP_to_smiles_and_names.pkl', 'rb'))
else:
    with open('/root/data/NP-MRD-dataset/NP-MRD_metadata/npmrd_natural_products.json') as f:
        data = json.load(f)
    # NP_to_inchi = dict([[d['accession'] , d['inchi']]for d in data['np_mrd']['natural_product']])
    NP_to_smiles = dict([[d['accession'] ,get_canonical_smiles(d)]for d in tqdm.tqdm(data['np_mrd']['natural_product'])])
    NP_to_names = dict([[d['accession'] ,d['name']]for d in data['np_mrd']['natural_product']])
    with open('/workspace/NP_to_smiles_and_names.pkl', 'wb') as file:
        pickle.dump([NP_to_smiles, NP_to_names], file)
    # exit(0)


#  get_nmr_tensors from txt file
def get_nmr_tensors(file_path):
    c_values = set()
    h_values = set()

    # print("opening file_path", file_path)
    with open(file_path, 'r') as file:
        # print('ys')
        for line in file:
            parts = line.strip().split(',')
            if len(parts) < 3 or not parts[2]:
                continue  # Skip lines that don't have enough parts or the third column is empty

            element, _, value, *_ = parts
            value = float(value)

            if element == 'C':
                c_values.add(value)
            elif element == 'H':
                h_values.add(value)

    # Convert sets to tensors
    c_tensor = torch.tensor(list(c_values), dtype=torch.float32)
    h_tensor = torch.tensor(list(h_values), dtype=torch.float32)
    # print("returning")
    return c_tensor, h_tensor
      

if __name__ == "__main__":
    """Used to add 1D NMR of molecules that we do NOT have 2D NMRs"""
    print("main")
    NP_MRD_FILES_dir = '/root/data/NP-MRD-dataset/NP-MRD-shift-assignments'
    npmrd_files_txt_only = os.listdir(NP_MRD_FILES_dir)
    
    oneD_only_name_and_smiles_mappings = {
        "train": [{},{}],
        "val": [{},{}],
        "test": [{},{}]
    }
    oneD_dataset_root_path = '/workspace/OneD_Only_Dataset'
    # print("start making directories")
    os.makedirs(oneD_dataset_root_path, exist_ok=True)
    os.makedirs(oneD_dataset_root_path+"/train/oneD_NMR", exist_ok=True)
    os.makedirs(oneD_dataset_root_path+"/val/oneD_NMR", exist_ok=True)
    os.makedirs(oneD_dataset_root_path+"/test/oneD_NMR", exist_ok=True)

            
    id_oneD_only_compound = 0
    collision = 0
    # for f in tqdm.tqdm(npmrd_files_txt_only):
    # print("start adding 1D NMR")
    for i, f in enumerate(tqdm.tqdm(sorted(npmrd_files_txt_only))):
        # if i % 100 == 0:  # Adjust the frequency based on your need
        #     snapshot = tracemalloc.take_snapshot()
        #     check_memory(snapshot)
        
        if f.split("_")[0]  not in NP_to_smiles:
            # something weird of NP-MRD that NP id(accession) not found in the json file
            continue
        
        smile =  NP_to_smiles[f.split("_")[0]]
        if smile: # sometimes it may just be None
            
            c_tensor, h_tensor = get_nmr_tensors(os.path.join(NP_MRD_FILES_dir, f))
            if smile not in canonical_smiles_to_path: # we do not have 2D NMR for this compound
                
                # 80% chance to send to training set
                rand_num =random.random()
                if  rand_num < 0.8:
                    split = "train" 
                elif rand_num < 0.9:
                    split = "val"
                else:
                    split = "test"
                    
                # also need to remember their chemical name and smiles
                name_mapping, smiles_mapping = oneD_only_name_and_smiles_mappings[split]
                name_mapping[id_oneD_only_compound] = NP_to_names[f.split("_")[0]]
                smiles_mapping[id_oneD_only_compound] = smile
                # print("will save at new split: ", split)
                # finally save  
                # print("save at", oneD_dataset_root_path+f"/{split}/oneD_NMR/{id_oneD_only_compound}.pt")
                torch.save([c_tensor, h_tensor], oneD_dataset_root_path+f"/{split}/oneD_NMR/{id_oneD_only_compound}.pt" )  
                id_oneD_only_compound+=1
            else: # we have 2D NMR for this compound
                
                file_path = canonical_smiles_to_path[smile]
                # print("save at", file_path)
                torch.save([c_tensor, h_tensor], file_path )
            # print("save finish")
           
        # if  i % 1000 == 0:  # Adjust based on how often you want to log memory usage
        #     current, peak = tracemalloc.get_traced_memory()
        #     logging.debug(f"Current memory usage: {current / 10**6}MB; Peak: {peak / 10**6}MB")
        
    # save the added chemical names and smiles
    for split in ["train", "val", "test"]:
        name_mapping, smiles_mapping = oneD_only_name_and_smiles_mappings[split]
        with open(f'/workspace/OneD_Only_Dataset/{split}/Chemical/index.pkl', 'w') as f:
            pickle.dump(name_mapping, f)
        with open(f'/workspace/OneD_Only_Dataset/{split}/SMILES/index.pkl', 'w') as f:
            pickle.dump(smiles_mapping, f)