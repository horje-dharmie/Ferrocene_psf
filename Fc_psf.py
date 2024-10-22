#!/home/damiojedeji/anaconda3/bin/python3
import numpy as np
import sys
import mdtraj as md
import random



####### Write PSF file so that bonds are display when using VMD ####################################

out_psf = open('FC.psf','w')

out_psf.write('PSF\n')
out_psf.write('       4 !NTITLE\n')
out_psf.write(' REMARKS Topology created\n')
out_psf.write(' REMARKS for Ferrocene \n\n')

#num_tween_mol = 2304
#num_toluene_mol=72000
num_ferrocene_mol= 1
#ncg_w= 104976
#ncgbeads = 29*num_tween_mol+3*num_toluene_mol+2*num_butanol_mol+ncg_w
# write atoms
num_ferrocene_atoms = num_ferrocene_mol * 23
out_psf.write('%8d !NATOM\n'%num_ferrocene_atoms)
mols=0
atoms=0

ferrocene_start=atoms
for imol in range(num_ferrocene_mol):
        mols+=1
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','X1', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','X2', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','CA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','HA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','CA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','HA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','CA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','HA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','CA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','HA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','CA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','HA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','CA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','HA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','CA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','HA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','CA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','HA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','CA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','HA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','CA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','HA', 0,0))) 
        atoms+=1
        out_psf.write('%8d    FC %-7d   FC%7s %-5s  %11.10f   %11.10f\n'%((atoms, mols, '   FC','Fe', 0,0))) 


out_psf.write('\n')
# write bonds
nbonds = num_ferrocene_mol*40

out_psf.write('%8d !NBOND: bonds\n'%nbonds)

bonds=0
for imol in range(num_ferrocene_mol):
    iat=ferrocene_start+imol*40+3
    bonds+=1
    out_psf.write('%8d%8d'%((iat,iat+1))) #3,4
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+2,iat+3)))#5,6
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+4,iat+5)))#7,8
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+6,iat+7)))#9,10
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+8,iat+9)))#11,12
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+10,iat+11)))#13,14
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+12,iat+13)))#15,16
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+14,iat+15)))#17,18
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+16,iat+17)))#19,20
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+18,iat+19)))#21,22
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat,iat+20)))#3,23
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+2,iat+20)))#5,23
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+4,iat+20)))#7,23
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+6,iat+20)))#9,23
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+8,iat+20)))#11,23
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+10,iat+20)))#13,23
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+12,iat+20)))#15,23
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+14,iat+20)))#17,23
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+16,iat+20)))#19,23
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+18,iat+20)))#21,23
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat,iat-2)))#3,1
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+2,iat-2)))#5,1
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+4,iat-2)))#7,1
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+6,iat-2)))#9,1
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+8,iat-2)))#11,1
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+10,iat-1)))#13,2
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+12,iat-1)))#15,2
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+14,iat-1)))#17,2
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+16,iat-1)))#19,2
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+18,iat-1)))#21,2
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat,iat+2)))#3,5
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+2,iat+4)))#5,7
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+4,iat+6)))#7,9
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+6,iat+8)))#9,11
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+8,iat)))#11,3
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+10,iat+12)))#13,15
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+12,iat+14)))#15,17
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+14,iat+16)))#17,19
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+16,iat+18)))#19,21
    if (bonds%4==0):
        out_psf.write('\n')
    bonds+=1
    out_psf.write('%8d%8d'%((iat+18,iat+10)))#21,13
    if (bonds%4==0):
        out_psf.write('\n')


out_psf.close()
