def get_s(chgcar='CHGCAR'):
    chg = VaspChargeDensity(chgcar)
    abc = np.shape(chg.chg[-1])
    dx = abc[0]/np.shape(chg.chg[0])[0]
    dy = abc[1]/np.shape(chg.chg[0])[1]
    dz = abc[2]/np.shape(chg.chg[0])[2]
    v_chg_gr = np.gradient(chg.chg[0],dx,dy,dz)
    chg_grad = np.sqrt( np.square(v_chg_gr[0]) + np.square(v_chg_gr[1]) + np.square(v_chg_gr[2]) )
    s = chg_grad / (2 * (3.* pi**2 * chg.chg[0])**1./3 * chg.chg[0] )
    return s

def backup_outcar():
    # it needs some improvement. Use with caution.
    import os
    import numpy as np
    from ase.io import read, write, Trajectory
    calc_id = int(np.loadtxt('db_id'))
    if os.path.isfile('OUTCAR') and not os.path.isfile(f'opt_id_{calc_id}.traj'):
        write(f'opt_id_{calc_id}.traj', read('OUTCAR',index=':'))
    elif os.path.isfile('OUTCAR') and os.path.isfile(f'opt_id_{calc_id}.traj'):
        t1 = Trajectory(f'opt_id_{calc_id}.traj', 'a')
        t2 = read('OUTCAR',index=':')
        for atoms in t2:
            t1.write(atoms)
        t1.close()
    return print('OUTCAR backup complete')

def plan_avg_ase(filename='LOCPOT'):
    import numpy as np
    from ase.calculators.vasp import VaspChargeDensity
    locpot = VaspChargeDensity(filename)
    # For LOCPOT files we multiply by the volume to get back to eV
    if 'LOCPOT' in filename:
        avg_c = [np.average(locpot.chg[0][:,:,i]*locpot.atoms[0].get_volume()) for i in range(0,np.shape(locpot.chg)[3]
)]
    else:
        avg_c = [np.average(locpot.chg[0][:,:,i]) for i in range(0,np.shape(locpot.chg)[3])]
    return avg_c
