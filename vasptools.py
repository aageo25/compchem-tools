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

def get_site_dcenter(vaspdos,site):
    """
    Fermi energy should be set to zero when defining vaspdos.
    Returs d-band center for spin up and down
    """
    from numpy import where,trapz
    
    dcenter = {'up': 0.0,'down': 0.0}
    for spin in dcenter.keys():
        # only energies below the fermi level
        e = vaspdos.energy[vaspdos.energy < 0]
        de = e[1] - e[0]
        # collect d states
        if spin == 'up':
            d = vaspdos.site_dos(site,8)[:len(e)]
            for i in range(10,18,2):
                d += vaspdos.site_dos(site,i)[:len(e)]
        elif spin == 'down':
            d = vaspdos.site_dos(site,9)[:len(e)]
            for i in range(11,19,2):
                d += vaspdos.site_dos(site,i)[:len(e)]
        # calculate d-band center
        num = trapz(d * e, dx=de)
        den = trapz(d, dx=de)
        dcenter[spin] = num / den
    return dcenter

def site_pdos(atoms,vaspdos):
    from ase.calculators.vasp import Vasp,VaspDos
    import numpy as np
    import matplotlib.pyplot as plt
    
    # collect the necessary information
    calc = Vasp(restart=True)
    atoms = calc.get_atoms()
    vaspdos = VaspDos(efermi=atoms.calc.fermi)
    
    colors = plt.get_cmap("tab10")
    
    # total dos
    energy = vaspdos.energy
    dos_up = vaspdos.dos[0]
    dos_down = vaspdos.dos[1]*(-1)
    
    # build the dictonary
    pdos = {}
    for atom in set(atoms.get_chemical_symbols()):
        pdos[atom] = {}
        for spin in ['up', 'down']:
            pdos[atom][spin] = np.zeros(energy.shape)
    
    # collect information from DOSCAR
    for symb in pdos.keys():
        for sp in pdos[symb].keys():
            for site in [atom.index for atom in atoms if atom.symbol == symb]:
                for orb in ["s", "py", 'pz', 'px', "dxy", "dyz", "dz2", "dxz", "dx2"]:
                    pdos[symb][sp] += vaspdos.site_dos(atom=site,orbital=f'{orb}-{sp}')
    
    # plot
    fig, ax = plt.subplots()
    # total
    ax.plot(energy, dos_up, label='Total', color='black')
    ax.plot(energy, dos_down, color='black')
    # element projected
    for index, symb in enumerate(pdos.keys()):
        for sp in pdos[symb].keys():
            if sp == 'up':
                ax.plot(energy, pdos[symb][sp], label=f'{symb}', color=colors(index))
                ax.fill_between(energy, pdos[symb][sp], color=colors(index), alpha=0.3)
            elif sp == 'down':
                ax.plot(energy, -pdos[symb][sp], color=colors(index))
                ax.fill_between(energy, -pdos[symb][sp], color=colors(index), alpha=0.3)
    
    # Fermi energy
    ax.axvline(0, ls='--', color='gray', label="Fermi energy")
    
    # Final adjustments
    ax.set_title('Element Projected DOS')
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('DOS (a.u.)')
    ax.legend()
    
    # return the graph
    return plt
