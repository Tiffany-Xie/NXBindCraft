###refold binder & target for evaluation
import numpy as np
import jax
import re
import pyrosetta as pr
from colabdesign import mk_af_model, clear_mem
from colabdesign.shared.utils import copy_dict
from functions.pyrosetta_utils import pr_relax, unaligned_rmsd, align_pdbs
from functions.biopython_utils import target_pdb_rmsd

# -------------------------------------------------------------------

# input: target pdb + binder sequence
#''' targete pdb, binder csv -> extract sequence -> '''
'''
complex pdb -> extract shorter one seq as binder -> extract longer one pdb as target -> 

'''

''' possible that fixing one structure during repredict limits realistic predictions (check)'''
mk_af_model.save_pdb
# ------------------------------------------------------------------- functions
main_folder = '/hpf/projects/mtyers/ningrui/NXBindCraft'
trial_name = 'testRefold'
def predict_complex(prediction_model, binder_sequence, complex_name, prediction_models, num_recycles_validation, filters, trial_name='testRefold'): # --> adapted from predict_binder_complex
    prediction_stats = {}

    # clean sequence
    binder_sequence = re.sub("[^A-Z]", "", binder_sequence.upper())

    # reset filtering conditionals
    pass_af2_filters = True
    # filter_failures = {}

    # start prediction per AF2 model, 2 are used by default due to masked templates
    for model_num in prediction_models:
        # check to make sure prediction does not exist already
        complex_pdb = os.path.join(main_folder, 'myTrials', trial_name, 'refoldPDB/refold', f'{complex_name}_model{model_num+1}.pdb')
        # if not os.path.exists(complex_pdb):
        prediction_model.predict(seq=binder_sequence, models=[model_num], num_recycles=num_recycles_validation, verbose=False)
        prediction_model.save_pdb(complex_pdb)
        prediction_metrics = copy_dict(prediction_model.aux['log'])
        
        # extract the statistics for the model
        stats = {
            'pLDDT': round(prediction_metrics['plddt'], 2), 
            'pTM': round(prediction_metrics['ptm'], 2), 
            'i_pTM': round(prediction_metrics['i_ptm'], 2), 
            'pAE': round(prediction_metrics['pae'], 2), 
            'i_pAE': round(prediction_metrics['i_pae'], 2)
        }
        

        # List of filter conditions and corresponding keys
        filter_conditions = [
            (f"{model_num+1}_pLDDT", 'plddt', '>='),
            (f"{model_num+1}_pTM", 'ptm', '>='),
            (f"{model_num+1}_i_pTM", 'i_ptm', '>='),
            (f"{model_num+1}_pAE", 'pae', '<='),
            (f"{model_num+1}_i_pAE", 'i_pae', '<='),
        ]

        for filter_name, metric_key, comparison in filter_conditions:
            threshold = filters.get(filter_name, {}).get("threshold")
            if threshold is not None:
                if comparison == '>=' and prediction_metrics[metric_key] < threshold:
                    pass_af2_filters = False
                    #filter_failures[filter_name] = filter_failures.get(filter_name, 0) + 1
                elif comparison == '<=' and prediction_metrics[metric_key] > threshold:
                    pass_af2_filters = False
                    #filter_failures[filter_name] = filter_failures.get(filter_name, 0) + 1

        #if not pass_af2_filters:
        #    break
        stats['pass_af2_filters'] = pass_af2_filters
        prediction_stats[model_num+1] = stats
        # end if not os.path.exists(complex_pdb)

    # skip update csv file with the failure counts
    
    # AF2 filters passed, contuing with relaxation
    for model_num in prediction_models:
        complex_pdb = os.path.join(main_folder, 'myTrials', trial_name, 'refoldPDB/refold', f'{complex_name}_model{model_num+1}.pdb')
        if prediction_stats[model_num+1]['pass_af2_filters']:
            complex_relaxed = os.path.join(main_folder, 'myTrials', trial_name, 'refoldPDB/refold_relaxed', f'{complex_name}_model{model_num+1}.pdb')
            pr_relax(complex_pdb, complex_relaxed)
        else:
            #if os.path.exists(complex_pdb):
            #    os.remove(complex_pdb)
            print('Complex prediction did not pass basic AF2 filters, but still keep it')
        
    return prediction_stats, pass_af2_filters

def predict_binder(prediction_model, binder_sequence, complex_name, gt_pdb, binder_chain, num_recycles_validation, prediction_models, trial_name='testRefold'):
    binder_stats = {}

    # prepare sequence for prediction
    binder_sequence = re.sub("[^A-Z]", "", binder_sequence.upper())
    prediction_model.set_seq(binder_sequence)

    for model_num in prediction_models:
        binder_alone_pdb = os.path.join(main_folder, "myTrials", trial_name, 'refoldPDB', 'refold_binder', f'{complex_name}_binder_model{model_num+1}.pdb')
        # if not os.path.exists(binder_alone_pdb):
        prediction_model.predict(models=[model_num], num_recycles=num_recycles_validation, verbose=False)
        prediction_model.save_pdb(binder_alone_pdb)
        prediction_metrics = copy_dict(prediction_model.aux['log'])

        # align binder model to trajectory binder
        align_pdbs(gt_pdb, binder_alone_pdb, binder_chain, 'A')

        # extract the statistics for the model
        stats = {
            'pLDDT': round(prediction_metrics['plddt'], 2), 
            'pTM': round(prediction_metrics['ptm'], 2), 
            'pAE': round(prediction_metrics['pae'], 2)
        }
        binder_stats[model_num+1] = stats
        # end if not os.path.exists(binder_alone_pdb)
    return binder_stats


# extract target pdb and binder sequence (longer = target; shorter = binder)
# TODO: whether the way separate binder/target affact prediction (i.e. whehter structure affact prediction e.g. # of helix)
from Bio import PDB
from Bio.PDB import PDBIO, Select
from Bio.SeqUtils import seq1
import os

class SelectChain(Select):
    def __init__(self, chain_id):
        super().__init__()
        self.chain_id = chain_id
    def accept_chain(self, chain):
        return chain.id == self.chain_id
    def accept_residue(self, residue):
        return PDB.is_aa(residue)
    

def separate_complex(complex_name):
    #complex_name = complex_pdb.removesuffix('.pdb')
    targets_path = '/hpf/projects/mtyers/ningrui/NXBindCraft/myTrials/testRefold/Targets'
    complex_pdb_path = os.path.join(targets_path, complex_name+'.pdb')

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('complex', complex_pdb_path)
    model = structure[0]

    chain_len = {}
    chain_seq = {}

    for chain in model:
        residues = [res for res in chain if PDB.is_aa(res)]
        seq = "".join(seq1(res.get_resname()) for res in residues)
        chain_id = chain.id
        chain_len[chain_id] = len(seq)
        chain_seq[chain_id] = seq

    binder_chain_id = min(chain_len, key=chain_len.get)
    target_chain_id = max(chain_len, key=chain_len.get)

    print(f'Target chain: {target_chain_id}')
    print(f'Binder chain: {binder_chain_id}')

    io = PDBIO()
    io.set_structure(structure)
    target_path = os.path.join(targets_path, f'{complex_name}_target.pdb')
    if os.path.exists(target_path):
        print(f'Target for complex {complex_name} already exists.')
    else:
        io.save(target_path, SelectChain(target_chain_id))

    return target_path, chain_seq[binder_chain_id], target_chain_id, binder_chain_id



# ------------------------------------------------------------------- input

target_seq = 'ETGASAKELACQEITVPLCKGIGYEYTYMPNQFNHDTQDEAGLEVHQFWPLVEIQCSPDLKFFLCSMYTPICLEDYKKPLPPCRSVCERAKAGCAPLMRQYGFAWPDRMRCDRLPEQGNPDTLCMDYERGTLEVLFQG'
binder_seq = 'MKILSKEEIEKRIRELERRLHELMTKLQPIFRDLPEEDRNEIIHLSLMNNDVWWDYRRNWETMSNEELQEKLEQLEEILKRFEEIVEKH'

params_path = '/hpf/projects/mtyers/ningrui//BindCraft'
target_pdb = '/hpf/projects/mtyers/ningrui/BindCraft/Targets/FZD8.pdb'
binder_csv = 'something'

# ------------------------------------------------------------------- target pdb + binder seq
num_recycles_validation = 3
multimer_validation = False
predicted_initial_guess = False
predict_bigbang = False

target_pdb, binder_seq, target_chain, binder_chain = separate_complex('5zng.pdb')

rm_template_seq_predict = False
rm_template_sc_predict = False


### compile model + prep inputs
# compile complex prediction model
complex_prediction_model = mk_af_model(protocol = 'binder',
                                       num_recycles = num_recycles_validation,
                                       data_dir = params_path,
                                       use_multimer = multimer_validation,
                                       use_initial_guess = predicted_initial_guess,
                                       use_initial_atom_pos = predict_bigbang)

# prepare input (complex)
complex_prediction_model.prep_inputs(pdb_filename=target_pdb, 
                                     chain='A', 
                                     binder_len=len(binder_seq),
                                     rm_target_seq=rm_template_seq_predict, 
                                     rm_target_sc=rm_template_sc_predict)


# compile binder monomer prediction model
binder_prediction_model = mk_af_model(protocol='hallucination',
                                      use_templates = False,
                                      initial_guess=False,
                                      use_initial_atom_pos=False,
                                      num_recycles=num_recycles_validation,
                                      data_dir=params_path,
                                      use_multimer=multimer_validation)

# prepare input (binder only)
binder_prediction_model.prep_inputs(length=len(binder_seq))

### predict
'''
HotspotRMSD -> unaligned RMSD of binder compared to original trajectory, in other words how far is binder in the repredicted complex from the original binding site
Target_RMSD -> RMSD of target predicted in context of the designed binder compared to input PDB
Binder_RMSD -> RMSD of binder predicted alone compared to original trajectory
'''
# complex predict
import json
filter_path = '/hpf/projects/mtyers/ningrui/NXBindCraft/settings_filters/default_filters.json'
with open(filter_path, 'r') as file:
    filters = json.load(file)

complex_statistics, pass_af2_filters = predict_complex(prediction_model=complex_prediction_model,
                                                       binder_sequence=binder_seq, 
                                                       complex_name='5ZNG',
                                                       prediction_models=[0,1],
                                                       num_recycles_validation=num_recycles_validation,
                                                       filters=filters)

# if not pass af2 filters, noted in dataframe; but also do scoring
# select the stats with highest plddt to continue scoring
pass_af2_keys = [k for k, v in complex_statistics.items() if v['pass_af2_filters']]
if pass_af2_keys:
    best_complex_pred_model = max(pass_af2_keys, key=lambda k: complex_statistics[k]['pLDDT'])
else:
    best_complex_pred_model = max(complex_statistics, key=lambda k: complex_statistics[k]['pLDDT'])

complex_name='5ZNG'
gt_complex_pdb = os.path.join(main_folder, 'myTrials/testRefold/Targets', '5zng.pdb')
# for model_num in prediction_models:
complex_pdb = os.path.join(main_folder, 'myTrials', trial_name, 'refoldPDB/refold', f'{complex_name}_model{best_complex_pred_model}.pdb')

pr.init('-ignore_unrecognized_res -ignore_zero_occupancy -mute all -holes:dalphaball /hpf/projects/mtyers/ningrui/NXBindCraft/functions/DAlphaBall.gcc -corrections::beta_nov16 true -relax:default_repeats 1')

if not os.path.exists(complex_pdb):
    print('Predicted complex structure not exists')
else:
    # unaligned RMSD to determine if binder is in the desinged binding site
    rmsd_site = unaligned_rmsd(gt_complex_pdb, complex_pdb, binder_chain, 'B')
    target_rmsd = target_pdb_rmsd(complex_pdb, target_pdb, target_chain)

rmsd_site
target_rmsd

### predict binder alone
# NOTE: paper only used template based model [0,1], here used [0-4]
binder_statistics = predict_binder(prediction_model=binder_prediction_model,
                                   binder_sequence=binder_seq,
                                   complex_name='5ZNG',
                                   gt_pdb=gt_complex_pdb,
                                   binder_chain=binder_chain,
                                   num_recycles_validation=num_recycles_validation,
                                   prediction_models=[0,1,2,3,4])
best_binder_pred_model = max(binder_statistics, key=lambda k: binder_statistics[k]['pLDDT'])
binder_pdb = os.path.join(main_folder, "myTrials", trial_name, 'refoldPDB', 'refold_binder', f'{complex_name}_binder_model{best_binder_pred_model}.pdb')
if not os.path.exists(binder_pdb):
    print('Predicted binder alone structure not exists')
else:
    rmsd_binder = unaligned_rmsd(gt_complex_pdb, binder_pdb, binder_chain, 'A')

rmsd_binder

# TODO: add another RMSD score to measure the difference between predicted complex binder and predicted binder alone
complex_pdb # predicted
binder_pdb # predicted
# align binder-binder first -> get score -> align back?a
align_pdbs(complex_pdb, binder_pdb, 'B', 'A')
rmsd_pre_binder = unaligned_rmsd(complex_pdb, binder_pdb, 'B', 'A')
rmsd_pre_binder
align_pdbs(gt_complex_pdb, binder_pdb, binder_chain, 'A')




# cn518




# ------------------------------------------------------------------- draft

import os
os.path.realpath('/hpf/projects/mtyers/ningrui/NXBindCraft/myScripts/draft.py')
os.getcwd()

# ------------------------------------------------------------------- target seq + binder seq

# compile prediction model

# af_multimer_model = mk_af_model(protocol = 'binder',
#                                 num_recycles = 3,
#                                 data_dir = params_path,
#                                 use_multimer = True,
#                                 use_initial_guess = False,
#                                 use_initial_atom_pos = False)
# 
# # prepare input
# # TODO: two cases 1) two pure sequence 2) target structure provided
# af_multimer_model.prep_inputs(pdb_filename=target_pdb, chain='A', binder_len=len(binder_seq),
#                               rm_target_seq=False, rm_target_sc=False)
# # complex_seq = f'{target_seq}:{binder_seq}'
# af_multimer_model.set_seq(complex_seq)
# 
# # predict
# af_multimer_model.predict()


