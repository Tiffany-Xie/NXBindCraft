## BindCraft Pipeline

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
        run MPNN to generate sequence for binders (default=20) vvv  ## why not use ralaxed trajectory??
            initialize mpnn model
            check whether keep interface from traj or rediign the whole binder by mpnn
            prepare inputs for mpnn
            sample mpnn sequence in parallel (20)
            return mpnn sequences
        create set of MPNN sequences with allowed aa composition and not already exists
        if mpnn sequences not empty
            increasing cycle if beta sheet traj
            compile complex prediction model
            prepare inputs (if...else...)
            compile binder monomer prediction model
            prepare inputs
            for mpnn seq in mpnn sequences
                add design to dict {seq, score, seqid}
                save fasta sequence if true
                predict mpnn redesigned binder *complex* using masked template target
                    for model in prediction models (2 total)
                        











