%%time
import mdtraj

new_traj = md.load_frame('lipids/sys.xtc', 0, top = 'lipids/sys.pdb')
new_traj = new_traj.remove_solvent()
new_traj = new_traj.atom_slice(new_traj.topology.select("resid 0"))
new_traj.save_pdb('lipid_iterload.pdb')
flag = 0
for chunk in md.iterload('lipids/sys.xtc', 1, top = 'lipids/sys.pdb'):
    chunk = chunk.remove_solvent()
    for lipid in range(chunk.n_residues):
        if lipid != 0 or flag != 0:
            buff = chunk.atom_slice(chunk.topology.select("resid %s" % lipid))
            new_traj = new_traj.join([buff])
    flag += 1
    if flag == 50:
        break
new_traj.superpose(new_traj, 0)
new_traj.save_xtc('2_128_iterload.xtc')