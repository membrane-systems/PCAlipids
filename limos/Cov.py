%%time
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

traj_frame = md.load_frame('one_lip_traj.xtc', 0, top = 'lipid_top.pdb')
Cov = [0.] * traj_frame.n_atoms * 3
Cov = [Cov] * traj_frame.n_atoms * 3
Cov = np.array(Cov)

mean_coor = traj_frame.xyz
mean_coor.shape
for frame in range(1, 500):
    traj_frame = md.load_frame('one_lip_traj.xtc', frame, top = 'lipid_top.pdb')
    mean_coor += traj_frame.xyz
mean_coor /= (frame + 1)

for frame in range(500):
    traj_frame = md.load_frame('one_lip_traj.xtc', frame, top = 'lipid_top.pdb')
    cent_coor = traj_frame.xyz - mean_coor
    for i in range(traj_frame.n_atoms):
        for j in range(traj_frame.n_atoms):
            Cov[3 * i][3 * j] += cent_coor[0][i][0] * cent_coor[0][j][0]
#             print(cent_coor[0][i][0] * cent_coor[0][j][0])
            Cov[3 * i][3 * j + 1] += cent_coor[0][i][0] * cent_coor[0][j][1]
            Cov[3 * i][3 * j + 2] += cent_coor[0][i][0] * cent_coor[0][j][2]
            
            Cov[3 * i + 1][3 * j] += cent_coor[0][i][1] * cent_coor[0][j][0]
            Cov[3 * i + 1][3 * j + 1] += cent_coor[0][i][1] * cent_coor[0][j][1]
            Cov[3 * i + 1][3 * j + 2] += cent_coor[0][i][1] * cent_coor[0][j][2]
            
            Cov[3 * i + 2][3 * j] += cent_coor[0][i][2] * cent_coor[0][j][0]
            Cov[3 * i + 2][3 * j + 1] += cent_coor[0][i][2] * cent_coor[0][j][1]
            Cov[3 * i + 2][3 * j + 2] += cent_coor[0][i][2] * cent_coor[0][j][2]
Cov /= (frame + 1)
print(np.linalg.eigvals(Cov))
plt.plot(np.linalg.eigvals(Cov)[:10])
plt.show()