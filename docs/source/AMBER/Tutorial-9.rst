.. note::

   Please report bugs, errors or enhancement requests through `Issue Tracker <https://github.com/SIRAHFF/documentation/issues>`_ or if you have a question about SIRAH open a `New Discussion <https://github.com/SIRAHFF/documentation/discussions>`_.
   
This tutorial shows how to use the SIRAH force field to perform a coarse grained (CG) simulation of a
Virus-like Particle in explicit solvent (called WatFour, WT4). The main references for
this tutorial are: `Darré et al. <https://pubs.acs.org/doi/abs/10.1021/ct100379f>`_, `Machado et al. <https://doi.org/10.1021/acs.jctc.9b00006>`__ and `Machado & Pantano  <https://academic.oup.com/bioinformatics/article/32/10/1568/1743152>`_. We strongly advise you to read these articles before starting the tutorial.

.. important::

    Check the :ref:`Setting up SIRAH <download amber>` section for download and set up details before starting this tutorial.
    Since this is **Tutorial 9**, remember to replace ``X.X`` and the files corresponding to this tutorial can be found in: ``sirah_[version].amber/tutorial/9/``


9.0. Structure preprocessing
____________________________

.. important::

    Check the :ref:`Setting up SIRAH <download amber>` section for download and set up details before starting this tutorial.
    Since this is **Tutorial 9**, remember to replace ``X.X`` and the files corresponding to this tutorial can be found in: ``sirah_[version].amber/tutorial/9/``

Create topology file and perform Energy Minimization for the capsomer:

.. code-block:: bash

  tleap -f 6OLA.leap
  cd run
  ./run_em.sh &  
  

Generate Capsid (complete VLP) using the capsomer previously generated:

.. code-block:: bash

  vmd -dispdev text -e capsid.tcl  
  

Some structures can bring with some ficticious "errors". Particularly, in cases such as PDB ID 6OLA which holds DNA fragments, each one containing a Phosphate atom from the next nucleotide, so this should be corrected erasing this atom.

.. code-block:: bash

  sed '/\(P\|OP1\|OP2\)[[:space:]]\+DC[[:space:]]\+[A-Z][[:space:]]\+1/d' 6ola-assembly.pdb > 6ola_edited.pdb


Fix pdb file correcting index numbers 0 to 99999 and include TER at the end of each chain, also add END at the end of the structure.

.. code-block:: bash

  vmd -dispdev text 6ola_edited.pdb -e fixpdb.tcl  


Assign protonation state with pdb2pqr using pH ~ 7.4

.. code-block:: bash

  python /home/msonora/Documents/apbs-pdb2pqr-master/pdb2pqr/pdb2pqr.py --with-ph=7.4 --ph-calc-method=propka --ff=amber --ffout=amber --chain --verbose 6ola_edited_OK.pdb 6ola_edited_OK.pqr  


9.1. Build CG representations
_____________________________

.. caution::

    The mapping to CG requires the correct protonation state of each residue at a given pH. We recommend using the `CHARMM-GUI server <https://www.charmm-gui.org/>`_ and use the **PDB Reader & Manipulator** to prepare your system. An account is required to access any of the CHARMM-GUI Input Generator modules, and it can take up to 24 hours to obtain one. 
    
    Other option is the `PDB2PQR server <https://server.poissonboltzmann.org/pdb2pqr>`_ and choosing the output naming scheme of AMBER for best compatibility. This server was utilized to generate the *PQR* file featured in this tutorial. Be aware that modified residues lacking parameters such as: MSE (seleno MET), TPO (phosphorylated TYR), SEP (phosphorylated SER) or others are deleted from the PQR file by the server. In that case, mutate the residues to their unmodified form before submitting the structure to the server.

Map the protonated atomistic structure of protein `1CRN <https://www.rcsb.org/structure/1CRN>`_ to its CG representation:   

.. code-block:: bash

  ./sirah.amber/tools/CGCONV/cgconv.pl -i 6ola_edited_OK.pdb -o 6ola_edited_OK_cg.pdb
  

The input file ``-i`` 1CRN.pqr contains the atomistic representation of `1CRN <https://www.rcsb.org/structure/1CRN>`_ structure at pH **7.0**, while the output ``-o`` 1CRN_cg.pdb is its SIRAH CG representation.

.. tip::

    This is the basic usage of the script **cgconv.pl**, you can learn other capabilities from its help by typing:

    .. code-block:: bash

        ./sirah.amber/tools/CGCONV/cgconv.pl -h 
        
.. note::

    **Pay attention to residue names when mapping structures from other atomistic force fields or experimental structures.** Although we provide compatibility for naming schemes in PDB, GMX, GROMOS, CHARMM and OPLS, there might always be some ambiguity in the residue naming, specially regarding protonation states, that may lead to a wrong mapping. For example, SIRAH Tools always maps the residue name “HIS” to a Histidine protonated at the epsilon nitrogen (:math:`N_{\epsilon}`) regardless the actual proton placement. Similarly, protonated Glutamic and Aspartic acid residues must be named “GLH” and “ASH”, otherwise they will be treated as negative charged residues. In addition, protonated and disulfide bonded Cysteines must be named “CYS” and “CYX” respectively. These kind of situations need to be carefully checked by the users. In all cases the residues preserve their identity when mapping and back-mapping the structures. Hence, the total charge of the protein should be the same at atomistic and SIRAH levels. You can check the following mapping file to be sure of the compatibility: ``sirah.amber/tools/CGCONV/maps/sirah_prot.map``.    

  
.. important::

    By default, charged termini are used, but it is possible to set them neutral by renaming the residues from **s**\[code\] to **a**\[code\] (Nt-acetylated) or **m**\[code\] (Ct-amidated) after mapping to CG, where \[code\] is the root residue name in SIRAH. For example, to set a neutral N-terminal Histidine protonated at epsilon nitrogen (:math:`N_{\epsilon}`) rename it from “sHe” to “aHe”.


Please check both PDB and PQR structures using VMD: 

.. code-block:: bash

  vmd -m 6ola_edited_OK.pqr 6ola_edited_OK_cg.pdb


Use tleap to calculate the charge of the VLP. Edit the 6OLA.leap file to match your requiriments

.. code-block:: bash

  tleap -f 6OLA.leap  
  
You can check the number of positive and negative atoms in the system:

.. code-block:: bash

  cat 6ola_edited_OK_cg.odb | grep -c "GC sL"
  cat 6ola_edited_OK_cg.odb | grep -c "GC sR"
  cat 6ola_edited_OK_cg.odb | grep -c "GC sE"
  cat 6ola_edited_OK_cg.odb | grep -c "GC sD"
  cat 6ola_edited_OK_cg.odb | grep -c "PX"
  

9.2. Wrapping up VLP system with SIRAH & packmol
_________________________________________________

Before packaging the VLP system is needed estimate the number of coarse grained water molecules (WT4) per compartment, lets say inner-virus and outer-virus.
To aproximate the number of solvent molecules in each compartment you can apply the a function that estimate the number of molecules in a given volume. Depending your system you could need to edit the radious in the layers_radius.dat file.

.. code-block:: bash

  awk -f calc_n.awk layers_radius.dat  
  
The outcome of this command is:

WT4_out = 20548
WT4_out = 6885

Naw use the software packmol to wrapp up the componente of system: VLP + solvent + ions.

.. seealso::

       The available electrolyte species in SIRAH force field are: ``Na⁺`` (NaW), ``K⁺`` (KW) and ``Cl⁻`` (ClW) which represent solvated ions in solution. One ion pair (e.g., NaW-ClW) each 34 WT4 molecules results in a salt concentration of ~0.15M (see :ref:`Appendix <Appendix>` for details). Counterions were added according to `Machado et al. <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953>`__.
       

In order to calculate exactly the number of solvent and ions molecules we do have to pass through packmol, we have to bare in mind the charge of the VLP. Particulartly, the 6OLA system charge is +420. To apply a charge balance we do need to add at least -420 negative ions in the inner core of the VLP. We are going to wrapp 620 ClW at the inner part of the virion to compensate VLP charge and lack of complete genome inside the viral particle. Entonces al final la capa externa de solvente tiene:

620 ClW
6885 - 620 = 6265 WT4

Now to reach 150 mM NaCl concentration on the outer layer of the WT4:

20548/34 = 604 ionic pairs

But remember the VLP has 420 positive charge. So 620 - 420 = 200
We can split this charge in order to do the system more balanced.
Entoces al final la capa externa de solvente tiene: 

604 + 100 = 704 NaW
604 - 100 = 504 ClW
20548 - 1208 = 19340 WT4

Edit the PCV2.pkm file to fix the number of solvent and ions according to your system.

.. code-block:: bash

  packmol < PCV2.pkm >> PCV2_packmol.log &  
  

From now on it is just normal Amber stuff!


5.2. Prepare LEaP input
_________________________

Use a text editor to create the file ``gensystem.leap`` including the following lines:

.. code-block:: console

    # Load SIRAH force field
    addPath ./sirah.amber
    source leaprc.sirah

    # Load model
    protein = loadpdb 1CRN_cg.pdb

    # Info on system charge
    charge protein  
    
    # Set S-S bridges
    bond protein.3.BSG protein.40.BSG
    bond protein.4.BSG protein.32.BSG
    bond protein.16.BSG protein.26.BSG

    # Add solvent, counterions and 0.15M NaCl
    # Tuned solute-solvent closeness for best hydration
    solvateOct protein WT4BOX 20 0.7
    addIonsRand protein NaW 22 ClW 22

    # Save Parms
    saveAmberParmNetcdf protein 1CRN_cg.prmtop 1CRN_cg.ncrst

    # EXIT
    quit

.. caution::

    Each disulfide bond must be defined explicitly in LEaP using the command bond, e.g.: “*bond unit.ri.BSG unit.rj.BSG*”. Where *ri* and *rj* correspond to the residue index in the topology file starting from 1, which may differ from the biological sequence in the PDB file. You can try the command *pdb4amber* to get those indexes from the atomistic structure, but be aware that it may not work if the Cysteine residues are too far away: 

    .. code-block:: bash

       pdb4amber -i sirah.amber/tutorial/5/1CRN.pqr -o 1CRN_aa.pdb && cat 1CRN_aa_sslink

    
.. seealso::

       The available electrolyte species in SIRAH force field are: ``Na⁺`` (NaW), ``K⁺`` (KW) and ``Cl⁻`` (ClW) which represent solvated ions in solution. One ion pair (e.g., NaW-ClW) each 34 WT4 molecules results in a salt concentration of ~0.15M (see :ref:`Appendix <Appendix>` for details). Counterions were added according to `Machado et al. <https://pubs.acs.org/doi/10.1021/acs.jctc.9b00953>`__.
       

5.3. Run LEaP 
____________________

Run the LEaP application to generate the molecular topology and initial coordinate files:

.. code-block:: bash

    tleap -f gensystem.leap

.. note::

    Warning messages about long, triangular or square bonds in ``leap.log`` file are fine and expected due to the CG topology of some residues.


This should create a topology file ``1CRN_cg.prmtop`` and a coordinate file ``1CRN_cg.ncrst``.

Use VMD to check how the CG model looks like and particularly the presence of disulfide bonds:

.. code-block:: bash

  vmd 1CRN_cg.prmtop 1CRN_cg.ncrst -e ./sirah.amber/tools/sirah_vmdtk.tcl


.. tip::

    VMD assigns default radius to unknown atom types, the script ``sirah_vmdtk.tcl`` sets the right
    ones, according to the CG representation. It also provides a kit of useful selection macros, coloring methods and backmapping utilities.
    Use the command ``sirah_help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.

5.4. Run the simulation
_______________________

Make a new folder for the run:

.. code-block:: bash

    mkdir -p run; cd run

The folder ``sirah.amber/tutorial/5/`` contains typical input files for energy minimization
(``em1_WT4.in`` and ``em2_WT4.in``), equilibration (``eq1_WT4.in`` and ``eq2_WT4.in``) and production (``md_WT4.in``) runs. Please check carefully the
input flags therein, in particular the definition of flag *chngmask=0* at *&ewald* section is **mandatory**.

.. tip::

    **Some commonly used flags in Amber**

   - ``-i``: Input file.
   - ``-o``: Output file.
   - ``-p``: Parameter/topology file.
   - ``-c``: Coordinate file.
   - ``-r``: Restart file.
   - ``-x``: Trajectory file.
   - ``-ref``: Reference file

.. caution::

    These input files are executed by the **GPU** implementation of ``pmemd.cuda``. Other available modules are ``sander`` or ``pmemd``, which are both **CPU** implementations of Amber.

.. note::

    The same input files can be used to run on CPU with the modules ``pmemd`` or ``sander``.
    
    
**Energy Minimization of side chains and solvent by restraining the backbone:**

.. code-block:: bash

    pmemd.cuda -O -i ../sirah.amber/tutorial/5/em1_WT4.in -p ../1CRN_cg.prmtop -c ../1CRN_cg.ncrst -ref ../1CRN_cg.ncrst -o 1CRN_cg_em1.out -r 1CRN_cg_em1.ncrst &
 
**Energy Minimization of whole system:**

.. code-block:: bash

    pmemd.cuda -O -i ../sirah.amber/tutorial/5/em2_WT4.in -p ../1CRN_cg.prmtop -c 1CRN_cg_em1.ncrst -o 1CRN_cg_em2.out -r 1CRN_cg_em2.ncrst &

**Solvent Equilibration (NPT):**

.. code-block:: bash

    pmemd.cuda -O -i ../sirah.amber/tutorial/5/eq1_WT4.in -p ../1CRN_cg.prmtop -c 1CRN_cg_em2.ncrst -ref 1CRN_cg_em2.ncrst -o 1CRN_cg_eq1.out -r 1CRN_cg_eq1.ncrst -x 1CRN_cg_eq1.nc &
  
.. caution::

    Option **restraintmask=:'1-46'** in input file ``eq1_WT4.in`` must be set specifically for each system to embrace all protein’s residues.

**Soft equilibration to improve side chain solvation (NPT):**

.. code-block:: bash

    pmemd.cuda -O -i ../sirah.amber/tutorial/5/eq2_WT4.in -p ../1CRN_cg.prmtop -c 1CRN_cg_eq1.ncrst -ref 1CRN_cg_eq1.ncrst -o 1CRN_cg_eq2.out -r 1CRN_cg_eq2.ncrst -x 1CRN_cg_eq2.nc &
  

**Production (1000ns):**

.. code-block:: bash

    pmemd.cuda -O -i ../sirah.amber/tutorial/5/md_WT4.in -p ../1CRN_cg.prmtop -c 1CRN_cg_eq2.ncrst -o 1CRN_cg_md.out -r 1CRN_cg_md.ncrst -x 1CRN_cg_md.nc &



5.5. Visualizing the simulation
________________________________

That’s it! Now you can analyze the trajectory.
Process the output trajectory to account for the Periodic Boundary Conditions (PBC):

.. code-block:: bash

      echo -e "autoimage\ngo\nquit\n" | cpptraj -p ../1CRN_cg.prmtop -y 1CRN_cg_md.nc -x 1CRN_cg_md_pbc.nc --interactive

Load the processed trajectory in VMD:

.. code-block::

    vmd ../1CRN_cg.prmtop ../1CRN_cg.ncrst 1CRN_cg_md.nc -e ../sirah.amber/tools/sirah_vmdtk.tcl

.. note::

     The file ``sirah_vmdtk.tcl`` is a Tcl script that is part of SIRAH Tools and contains the macros to properly visualize the coarse-grained structures in VMD. Use the command ``sirah-help`` in the Tcl/Tk console of VMD to access the manual pages. To learn about SIRAH Tools' capabilities, you can also go to the :ref:`SIRAH Tools tutorial <SIRAH tools>`.
