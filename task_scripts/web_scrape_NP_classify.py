import  requests, json,pickle,tqdm, os

def get_superclass_and_glycoside(smiles):
    # smiles = 'CC1C(O)CC2C1C(OC1OC(COC(C)=O)C(O)C(O)C1O)OC=C2C(O)=O'
    url = f"https://npclassifier.gnps2.org/classify?smiles={smiles}"
    response = requests.get(url)
    json_dat = json.loads(response.content)
    print(json_dat)
    superclass_results = json_dat['get_superclass_and_glycoside'][0]
    isglycoside = json_dat['isglycoside']
    return superclass_results, isglycoside


# for split in ['val', 'test', "train"]:
#     superclass_results = {}
#     smiles_pkl = pickle.load(open(f'/workspace/SMILES_dataset/{split}/SMILES/index.pkl','rb'))
#     for k,v in tqdm.tqdm(smiles_pkl.items()):
        
#         superclass_results[k] = get_superclass_and_glycoside(v)
#     os.makedirs(f'/root/gurusmart/MorganFP_prediction/reproduce_previous_works/reproducing_deepsat/{split}/Superclass',exist_ok=True)            
#     pickle.dump(superclass_results,open(f'/root/gurusmart/MorganFP_prediction/reproduce_previous_works/reproducing_deepsat/{split}/Superclass/index.pkl','wb'))
        
print(get_superclass_and_glycoside('CC1C(O)CC2C1C(OC1OC(COC(C)=O)C(O)C(O)C1O)OC=C2C(O)=O'))