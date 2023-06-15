import numpy as np, subprocess
import sys
if (len(sys.argv)!=6):
    print("syntax:  material_budget.py input.xml eta")
    exit()
xml_file=sys.argv[1]
eta = float(sys.argv[2])
phimin=float(sys.argv[3])
phimax=float(sys.argv[4])
nphi=int(sys.argv[5])

useTilt=False
if useTilt:
    tilt=-.025
else:
    tilt = 0

tag=xml_file.split("/")[-1].split(".")[0]
    
#prints the nuclear interaction lengths and radiation lengths for a given eta as a function of phi
for theta in [2*np.arctan(np.exp(-eta))]:
    phis = np.linspace(phimin, phimax, nphi)
    nils = []
    rads=[]
    for phi in phis:
        Z=500 if eta>0 else -315
        uxp = np.sin(theta)*np.cos(phi)
        uyp = np.sin(theta)*np.sin(phi)
        uzp = np.cos(theta)
        ux = uxp*np.cos(tilt)+uzp*np.sin(tilt)
        uy = uyp
        uz = uzp*np.cos(tilt)-uxp*np.sin(tilt)
        command=f"print_materials {xml_file} 0 0 0 {Z*ux/uz} {Z*uy/uz} {Z}"
        #print(command)
        #command += " | grep 'Integrated interaction lengths' | awk '{print $5;}'"
        result = subprocess.run(command.split(), stdout=subprocess.PIPE, text=True, stderr=subprocess.DEVNULL)
        nil=-1
        rad=-1
        for line in result.stdout.split("\n"):
            if "Integrated interaction lengths" in line:
                nil=float((line.strip().split()[4]))
                #break;
            if "Integrated radiation lengths"  in line:
                rad=float((line.strip().split()[4]))
                #break;
        nils.append(nil)
        rads.append(rad)
    print(f"nils[('{tag}',{eta})] = ",nils)
    print(f"rads[('{tag}',{eta})] = ",rads)
        #val = float(result.stdout)
        #print(theta, phi, val)
