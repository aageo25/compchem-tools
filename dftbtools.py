def add_vasp_structure(dbname, folder, system_id=''):
   """Add a VASP structure into a database.

   Keyword arguments:
   dbname -- database filename
   folder -- path with VASP outputs
   system_id -- type of calculation you performed

   """
   import sys
   import os
   import re
   import ase.db
   from ase.io import read
   from ase.io.vasp import read_vasp_out


   start_dir = os.path.realpath(".")
   db = ase.db.connect(dbname)
   os.chdir(folder)

   structure = read_vasp_out("OUTCAR")

   fin = open("INCAR"); INCAR = fin.readlines(); fin.close()
   fin = open("KPOINTS"); KPOINTS = fin.readlines(); fin.close()
   fin = open("POTCAR");  POTCARf = fin.read();      fin.close();
   POTCAR = re.findall('TITEL  = (.*)',POTCARf)
   fin = open("CONTCAR"); CONTCAR= fin.readlines(); fin.close()

   datain = {"INCAR" : INCAR,
             "KPOINTS" : KPOINTS,
             "POTCAR_TITEL" : POTCAR,
             "CONTCAR" : CONTCAR}

   str_path = os.path.realpath(".")
   #str_system = str(type_system)
   os.chdir(start_dir)

   row_id = db.write(structure,data=datain, structure_path=str_path,system_id=system_id)
   print ("Structure added with id:", row_id)
   return row_id


def db_read_row(dbname, id):
    """Reads a specific row in a ASE database.

    Keyword arguments:
    dbname -- database name.
    id -- id number of the row to read.

    """
    import ase.db

    db = ase.db.connect(dbname)

    # Read all information for one row in the database
    row = db.get(id=id)

    structure = row.toatoms()
    structure.row = row
    return structure

def have_vasprun(location):
    import os
    """
    Checks if the vasprun exists.
    """

    file_list= os.listdir(location)
    flag = False

    if "vasprun.xml" in file_list:    #checks to see if vasprun.xml is in the list of files
        flag = True
        return flag
    else:
        return flag

def write_geo(atoms_obj, identifier, file_out='atoms_out.geo'):
    """Create a geo file from an ASE atoms object

    Keyword arguments:
    atoms_obj -- One ASE atoms object
    identifier -- unique id to parse the .geo file and the trainset.in
    output -- name of the geo file
    """
    struct = atoms_obj.copy()
    identifier = str(identifier)
    #file_out = output+'.geo'
    fout = open(file_out,"a")
    fout.write("XTLGRF 200\n")
    fout.write("DESCRP "+identifier+"\n")
    fout.write('CRYSTX ')
    for i in struct.get_cell_lengths_and_angles():
        fout.write(str(i)+' ')
    fout.write('\n')
    for i in range(len(struct.get_chemical_symbols())):
        fout.write(str('HETATM  '+str(i+1)+' '+str(struct.get_chemical_symbols()[i])+'     '+str(struct.positions[i])[1
:-1])+'\n')
    fout.write('END'+2*'\n')
    fout.close()
    return print("ID", identifier,"added in",file_out)

def write_trainset(atoms_obj, ref_obj, atoms_id, ref_id, weight='1.0', file_out='trainset.in'):
    """ Write/update the ENERGY block in the trainset.in file for reaxff

    Keyword arguments:
    atoms_obj -- list of atoms objects to include in the file
    ref_obj -- single atoms object to use as a reference for the dataset
    atoms_id -- identifier in the geo file for each atoms object. Must follow indexing in the atom_obj list
    ref_id -- identifier in the geo file for the reference atoms object.
    weight -- weight for the dataset

    """
    try:
        readfile = open(file_out,'r')
        info = 'updating'
        info_f = 'updated'
    except:
        readfile = open(file_out,'w')
        info = 'creating'
        info_f = 'created'
        readfile.writelines(['ENERGY\n','ENDENERGY\n'])
        readfile.close()
        readfile = open(file_out,'r')
    print(info,file_out)
    lines = readfile.readlines()
    lines = lines[1:-1]
    readfile.close()

    fout = open(file_out, "w")
    fout.write("ENERGY\n")
    fout.writelines([item for item in lines])
    for i in range(len(atoms_obj)):
        DeltaE = atoms_obj[i].get_total_energy() - ref_obj.get_total_energy()
        if str(atoms_id[i]) != str(ref_id):
            fout.write(str(weight+'    +  '+str(atoms_id[i])+'  -  '+str(ref_id)+'    '+str(DeltaE))+'\n')
    fout.write("ENDENERGY\n")
    fout.close()
    return print("File", file_out,"succesfully",info_f)

def prep_runfolders(dbname,query):
    import os
    import shutil
    from ase.db import connect
    
    db = connect(dbname)
    
    prevdir = os.getcwd()
    
    for row in db.select(query):
        dir = str(row.id)
        try:
            os.mkdir(dir)
        except FileExistsError:
            print(f'Keeping folder {dir}')
        else:
            print(f'Creating folder {dir}')
        os.chdir(dir)
        try:
            os.symlink('../run.sh', 'run.sh')
        except:
            pass
        with open('db_id', 'w') as out:
            out.write(dir)
        os.chdir(prevdir)
    return print('Done')
