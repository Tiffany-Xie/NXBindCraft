###refold binder & target for evaluation
import numpy as np
import jax
from colabdesign import mk_af_model, clear_mem
# -------------------------------------------------------------------

# input: target pdb + binder sequence
''' targete pdb, binder csv -> extract sequence -> '''
''' possible that fixing one structure during repredict limits realistic predictions (check)'''

# ------------------------------------------------------------------- input

target_seq = 'ABCD'
binder_seq = 'OIUY'
params_path = '/hpf/projects/mtyers/ningrui//BindCraft/params/'
target_pdb = 'something'
binder_csv = 'something'

# -------------------------------------------------------------------

# compile prediction model

af_multimer_model = mk_af_model(protocol = 'binder',
                                num_recycles = 3,
                                data_dir = params_path,
                                use_multimer = True,
                                use_initial_guess = False,
                                use_initial_atom_pos = False)

# prepare input
# TODO: two cases 1) two pure sequence 2) target structure provided
# af_multimer_model.prep_inputs(pdb_filename=1) = 
complex_seq = f'{target_seq}:{binder_seq}'
af_multimer_model.set_seq(complex_seq)


# predict

