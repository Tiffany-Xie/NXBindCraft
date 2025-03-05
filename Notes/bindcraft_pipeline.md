## BindCraft Pipeline

'''
seed, length, helicity value
trajectory (binder hallucination)
    initialize binder hallucination model
    prepare input for binder design
    update weights put on different metric
    check if add additional losses
    cahculate # mutations based on binder length
    Starting Design (2,3 stages, mcmc, greedy...) [4 stages]
        Testing Logits vvv (I)
        design logits for 50 iterations
        get initial plddt (after 50 iters) 
        if initial plddt > 0.65
            if beta sheet traj detected
                increase number of recycles
            if logit iterations left
                design logits the rest iterations
                update logit plddt
            Softmax Optimizations vvv (II)
            design soft
            get softmax plddt
            if softmax plddt > 0.65 
                One-hot Optimisation vvv (III)
                design hard
                get onehot plddt
                if onehot plddt > 0.65
                    PSSM Semigreddy Optimisation vvv (IV)
                    design pssm semigreddy
                else: onehot traj plddt too low
            else: softmax plddt too low
        else: initial plddt too low
    else: Error, no valid design model
    get final plddt; save model pdb; traj terminate = ''
    
    get C alpha clashes
    if ca clashes > 0
        terminate = 'clashing'
    else:
        if final plddt < 0.7
            terminate = 'lowConfidence'
        else:
            get binder contacts (identify residues at the binder interface within a range)
            get # binder contacts
            if # binder contacts < 3
                terminate = 'lowConfidence'
            else:
                terminate = ''
    if terminate != ''
        move traj pdb to trajectory/terminate path
    get sampled sequence for plotting
    save design trajectory plots
    save design animations
    save trajectory pickle (default=false)
    return af_model

get trajectory metrics (olddt, ptm, iptm, pae, i_pae)
if terminate = ''
    relaxed trajectory
    get clashes before relax
    get clashes after relax
    get secondary structure content of starting trajectory binder and interface
    get interface score, interface AA, interface residue for relaxed trajectory
    get binder sequence (best)
    analyze sequence get notes
    get RMSD for input PDB and trajectory target
    summarize trajectory statistics
    save trajectory statistics to CSV
    
    Start MPNN vvv
        run MPNN to generate sequence for binders (default=20) vvv  ## why not use ralaxed trajectory?? A: because training using cryst structure
            initialize mpnn model
            check whether keep interface from traj or rediign the whole binder by mpnn
            prepare inputs for mpnn
            sample mpnn sequence in parallel (20)
            return mpnn sequences
        create set of MPNN sequences with allowed aa composition and not already exists
        if mpnn sequences not empty
            increasing cycle if beta sheet traj
            compile complex prediction model ('binder', #recycle, afparams, x multimrValid,...)
            if use initial guess
                prepare input (trap pdb, 'A', 'B', binderlen, useBinderTemplate, rm target template for repred, rm sidechain from target for repred...)
                # Why need binder length if this is pure prediction ???
            else
                prepare input (starting pdb, 'B', binderlen, rm target template for repred, rm sidechain from target for repred)
            compile binder monomer prediction model ('hallucination', #recycles, af params, x multimrValid,...)
            prepare inputs (binder length)
            for mpnn seq in mpnn sequences
                add design to dict {seq, score, seqid}
                save fasta sequence if true
                predict mpnn redesigned binder *complex* using masked template target
                    for model in prediction models (2 total)
                        predict complex using prediction model
                        get prediction metric (plddt, ptm, i_ptm, pae, i_pae)
                        extract the stats of the model (plddt, ptm, i_ptm, pae, i_pae)
                        add stats to prediction stats dict (2 total)
                        list the initial filter conditions
                        for filter in filter conditions
                            check if the prediction stats pass the filter
                        if not pass filter: break
                    update failure csv
                    for mudel in prediction model (2 total) # relaxation
                        relax complex
                    return complex prediction stats dict and T/F pass af2 initial filters
                if not pass af2 initial filters: continue
                for model in prediction models (2 tot) # calculate stats
                    get clashes before and after relax
                    get mpnn interface scores, mpnn interface AA and mpnn interface residues fore relaxed mpnn complex
                    get secondary structure content of mpnn binder
                    get unaligned RMSD of traj binder vs. mpnn binder  # check if binder at the binding site
                    get RMAS of mpnn target vs. input PDB
                    add additional stats to mpnn complex predicion stats  # long list; refer to code
                calculate complex averages
                predict mpnn *binder alone* in single sequence mode
                    for model in prediction model
                        predict binder alone
                        get prediction metric (plddt, ptm, pae)
                        align predicted binder with traj binder
                        extract stats of the model (plddt, ptm, pae)
                        add stats to prediction binder stats dict
                    return binder stats dict
                for model in prediction model
                    get RMSD of mpnn predict binder vs. traj binder
                    add binder RMSD to binder stats dict
                calculate binder averages
                analyze mpnn sequence to make sure no cysteins and contains residues that absorb UV (note)
                insert mpnn design stats to CSV (non stats data + mpnn_complex + binder_alone)
                find the best model by plddt
                run design data against filter thresholds
                if filter conditions = True
                    good...
                else:
                    update failure csv
                if enought mpnn seq of the same traj: break
            if # accept mpnn >= 1: print mpnn designs passing filters
            else: no accepted mpnn design in this traj
        else: duplicate mpnn designs sampled with different traj, skipping current traj optimization
        remove unrelaxed designs traj PBD
        measure time
    analyse the rejetion rate of traj
    if acceptance < acceptance rate: break
trajectory + 1
'''
### Yeah the end!!!

                    
                    






