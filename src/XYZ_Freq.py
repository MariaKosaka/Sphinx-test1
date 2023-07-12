"""Frequency calculation module

This program can analyse vibration modes of molecule, visualize 
the normal mode coordinates, and calculate IR spectrum. 
Before executing this program, please prepare 
the chemical structure files in the current directory: 
<jobname>.xyz is necessary. 
If you want to visualize normal mode vectors, 
<jobname>.vesta is also required.

Usage:
    python XYZ_Freq.py

    >>> Enter the job name: <jobname>

"""


from pyscf import gto, scf
import numpy as np
import pymatgen.core.periodic_table as mg
import pylab as plt

class Hessian:
    """Hessian calculation class.
    
    Attributes:
        coord (numpy.ndarray): An array of atom coordinates [Angstroms].
        atmlist (list): A list of atoms (String).

        These are created from <jobname>.xyz via xyzImp function.
        When atmlist = ["C", "O"], 
        coord = [C_x (x coordinate for C atom), C_y, C_z, O_x, O_y, O_z].

        strc (String): A string defined with coord and atmlist.

    Usage:
        h = Hessian(coord, atmlist)
        hessian = h.hessian()


    """



    #initial var
    _d = 1.0e-5

    #Define Functions
    def __init__(self, coord: np.ndarray, atmlist: list):
        self.coord = coord
        self.dim = len(coord)
        self.atmlist = atmlist
    
        strc = ''
        for i in range(len(atmlist)):
            strc = strc + f'{atmlist[i]} {coord[3*i]} {coord[3*i+1]} {coord[3*i+2]}; '

        if strc[-2] == ';':
            strc = strc[:-2]

        self.strc = strc
     
    def energy(self) -> float:
        """Calculate SCF energy [Hartree].
        
        Returns:
            float: the SCF energy of the molecule

        """
        mol = gto.M(atom=self.strc, basis='cc-pvdz', verbose=0)
        mf = scf.RHF(mol)
        ene = mf.kernel()
        return ene

    def forceCalc(self, n: int) -> float:
        """Numerical derivative of Energy with nth coordinate. (dE/dx_n)
        
        Args:
            n (Int): the suffix of derivative variable

        Returns:
            float: the derivative of molecule's energy

        """
        coord = self.coord
        atmlist = self.atmlist
        dim = self.dim
        diff = np.zeros(dim, dtype=float)
        diff[n] = self._d
        __large = Hessian(coord + diff, atmlist)
        __small = Hessian(coord - diff, atmlist)
        force = (__large.energy()-__small.energy())/(2*self._d)
        return force
    
    def secDrv(self, n: int, m: int) -> float:
        """Second derivative of Energy with n,mth coordinates. (d^2 E/dx_n dx_m)

        Args:
            n, m (Int): the suffix of derivative variables

        Returns:
            float: the second derivative of molecule's energy

        """
        coord = self.coord
        atmlist = self.atmlist
        dim = self.dim
        diff = np.zeros(dim, dtype=float)
        diff[m] = self._d
        __large = Hessian(coord + diff, atmlist)
        __small = Hessian(coord - diff, atmlist)
        secDeriv = (__large.forceCalc(n)-__small.forceCalc(n))/(2*self._d)
        return secDeriv
    
    def hessian(self) -> np.matrix:
        """Calculate Hessian

        Returns:
            numpy.matrix: the Hessian of the matrix

        """
        coord = self.coord
        hess = np.empty((self.dim,self.dim), dtype=float)
        for i in range(self.dim):
            for j in range(i + 1):
                hess[i, j] = self.secDrv(i, j)
                hess[j, i] = hess[i, j]
        return hess


class Harmonic:
    """Frequency analysis class.
    
    Attributes:
        hess (numpy): Hessian, generated via Hessian.hessian()
        atmlist (list): A list of atoms (String), the same as Hessian.

        eig (list): Eigenvalues and eigenvectors of mass-weighted Hessian, 
                    calculated via Harmonic.frequencies() function.

    Usage:
        harmo = Harmonic(hess, atmlist)
        print(harmo)
        frequencies = harmo.eig[0]
        normalModeVectors = harmo.eig[1]

    """

    #constants
    _amu = 1.66e-27 #amu -> kg
    _hart = 4.36e-18 #Hartree -> J
    _pi = np.pi
    _c = 3.0e+10 #cm/s
    _angs = 1e-10 #Angstrom -> Meter

    def __init__(self, hess: np.matrix, atmlist:list):
        self.hess = hess
        self.atmlist = atmlist
        self.dim = len(atmlist) * 3
        self.eig = self.frequencies()

    def massWeighted(self) -> np.matrix:
        """Convert cartesan-based Hessian to Mass-Weighted Hessian.

        Returns:
            numpy.matrix: mass-weighted Hessian

        """
        mass = []
        for i in range(self.dim):
            atm = mg.Element(self.atmlist[i//3])
            atmMass = atm.atomic_mass * self._amu
            mass.append(atmMass)

        hess = self.hess * self._hart / (self._angs**2)
        mwHess = hess
        for i in range(self.dim):
            for j in range(i + 1):
                mwHess[i, j] = hess[i, j] / ((mass[i]*mass[j])**0.5)
                mwHess[j, i] = mwHess[i, j]
        return mwHess

    def frequencies(self) -> list:
        """Calculate eigenvalues and eigenvectors for mass-weighted Hessian.

        Returns:
            list: Eigenvalues and eigenvectors of mass-weighted Hessian

"""
        eig = np.linalg.eig(self.massWeighted())
        return eig

    def __str__(self):
        output = "Freq. Calculation Processing ...\n"
        freq = self.eig
        nmode = len(freq[0])
        for i in range(nmode):
            omega2 = freq[0][i]
            norMode = freq[1][:, i]

            notion = "*******************************"
            if omega2 < 0:
                notion = "*Warning: Imaginary Frequency!*"

            frequency = (abs(omega2)**0.5) / (2*self._pi*self._c) 

            output = output + f"\nMODE #{i+1}  {notion}\n"
            output = output + f"|Frequency| [cm**-1] = {frequency: 7.4f}\n"
            output = output + f"Normal Mode Coordinates [mwc]\n"

            for i in range(len(self.atmlist)):
                output = output + f"{self.atmlist[i]}  {norMode[3*i]: 7.4f}  {norMode[3*i+1]: 7.4f}  {norMode[3*i+2]: 7.4f}\n"

            output = output + "****************************************\n\n"
                       
        output = output + f"Freq. Calculation successfully terminated.\n\n"
        output = output + f"****************************************\n"
        return output


class IR(Hessian, Harmonic):
    """IR intensity calculation class.
    
    Attributes:
        eig (list): Eigenvalues and eigenvectors of mass-weighted Hessian, 
                    calculated via Harmonic.frequencies() function.
        coord (numpy.ndarray): An array of atom coordinates.
        atmlist (list): A list of atoms (String), the same as Hessian.
        intensity (numpy.ndarray): An array of IR intensity,
                  calculated via IR.xyzIntensity() or IR.normodeIntensity()

    Usage:
        ir = IR(eig, coord, atmlist)
        intensity = ir.xyzIntensity()
        print(ir)

    Note:
        Execute IR.xyzIntensity() or IR.normodeIntensity() before print(IR)!


    """

    def __init__(self, eig: np.ndarray, coord: np.ndarray, atmlist:list):
        self.freq = eig
        self.frequencies = np.abs(self.freq[0])**0.5 / (2*Harmonic._pi*Harmonic._c)
        Hessian.__init__(self, coord, atmlist)
        self.dim = len(atmlist)*3

    def dipoleMoment(self) -> np.ndarray:
        """Calculate Dipole Moment

        Returns:
            numpy.ndarray: dipolemoment [x, y, z]

        """
        mol = gto.M(atom=self.strc, basis='cc-pvdz', verbose=0)
        mf = scf.RHF(mol)
        mf.kernel()
        dipole = mf.dip_moment()
        return dipole

    def xyzIntensity(self) -> np.ndarray:
        """Calculate IR intensity of each modes,
           with unitary transformation of xyz derivative of dipole moment.

        Returns:
            numpy.ndarray: an array of IR intensity of each modes. 

        """
        coord = self.coord
        atmlist = self.atmlist
        dim = self.dim
        xyzIntensity = []
        for n in range(dim):
            diff = np.zeros(dim, dtype=float)
            diff[n] = self._d
            __large = IR(self.freq, coord + diff, atmlist)
            __small = IR(self.freq, coord - diff, atmlist)
            grad = (__large.dipoleMoment()-__small.dipoleMoment())/(2*self._d)
            peak = np.linalg.norm(grad, ord=2)**2
            xyzIntensity.append(peak)
        xyzIntensity = np.array(xyzIntensity)
        unitary = self.freq[1]
        invUni = np.linalg.inv(unitary)
        transformedIntensity = np.abs(np.dot(unitary, xyzIntensity)) / np.abs(self.frequencies) 
        #transformedIntensity = np.abs(np.dot(invUni, np.dot(xyzIntensity,unitary))) / np.abs(self.frequencies) 
        self.intensity = transformedIntensity
        return transformedIntensity

    def normodeIntensity(self) -> np.ndarray:
        """Calculate IR intensity of each modes,
           with derivative of dipole moment with normal mode coordinate.

        Returns:
            numpy.ndarray: an array of IR intensity of each modes. 

        """
        coord = self.coord
        atmlist = self.atmlist
        dim = self.dim

        intensity = []
        for n in range(len(self.freq[0])):
            frequency = self.frequencies[n]
    
            normodes = self.freq[1]
            nmode = normodes[:, n]
            diff = self._d*nmode/np.linalg.norm(nmode)
            __large = IR(self.freq, coord + diff, atmlist)
            __small = IR(self.freq, coord - diff, atmlist)
            grad = (__large.dipoleMoment()-__small.dipoleMoment())/(2*self._d)
            peak = np.linalg.norm(grad, ord=2)**2 / abs(frequency)
            intensity.append(peak)
        intensity = np.array(intensity)
        self.intensity = intensity
        return intensity

    def __str__(self):
        intensities = self.intensity
        frequencies = self.frequencies 
        output = "\nIR Spectrums\n#mode Frequency[cm**-1]  IR intensity\n"
        for mode in range(len(frequencies)):
            output = output + f"{mode+1: >5d} {frequencies[mode]: 17.8f}  {intensities[mode]: 12.6f}\n"
        return output


class Plot:
    """Plot and Visualization class.

    Attributes: 
        eig (list): Eigenvalues and eigenvectors of mass-weighted Hessian, 
                    calculated via Harmonic.frequencies() function.
        job (String): the name of molecule you want to calculate.
        intensity (numpy.ndarray): An array of IR intensity,
                  calculated via IR.xyzIntensity() or IR.normodeIntensity()
        atmlist (list): A list of atoms (String), the same as Hessian.

    Usage:
        Plot(eig, job, intensity, atmlist)
        Then IR spectrum graph
         and VESTA files for each normal mode vector
         which visualize them on the molecule are generated. 
  
    """

    _c = 3.0e+10 #cm/s
    _pi = np.pi

    def __init__(self, eig:np.ndarray, job:str, intensity:np.ndarray, atmlist:list):
        self.freq = eig
        self.job = job
        self.natm = len(atmlist)
        self.intensity = intensity 
        self.vecviz(1)
        self.irSpectrum()

    def vecviz(self, scale:float):
        """Visualize the Normal Mode Vectors on VESTA file.

        It generates <jobname>_vector_mode<number>.vesta file
        for each oscillation mode.
        
        Args:
            scale (float): The length scale of the vectors

        """
        try:
            normodes = self.freq[1]
            
            for mode in range(len(normodes)):
                nmode = normodes[:, mode]
                vectr = "VECTR\n"
                vectt = "VECTT\n"
                for atm in range(self.natm):
                    vectr = vectr + f"   {atm+1}    {nmode[3*atm]*scale}    {nmode[3*atm+1]*scale}    {nmode[3*atm+2]*scale} 0\n"
                    vectr = vectr + f"    {atm+1}   0    0    0    0\n"
                    vectr = vectr + f" 0 0 0 0 0\n"

                    vectt = vectt + f"   {atm+1}  0.300 255   0   0 2\n"
                
                vectr = vectr + " 0 0 0 0 0\n"
                vectt = vectt + " 0 0 0 0 0\n"
        
                replace = vectr + vectt
                original = "VECTR\n"
                original = original + f" 0 0 0 0 0\n"
                original = original + f"VECTT\n"
                original = original + f" 0 0 0 0 0\n"
    
                with open(self.job+'.vesta', 'r') as f:
                    data_lines = f.read()
                data_lines = data_lines.replace(original, replace)
                with open(self.job+'_vector_mode'+str(mode+1)+'.vesta', mode="w") as f:
                    f.write(data_lines)

        except FileNotFoundError:
            print("There's no " + self.job +".vesta file. Please check the file name.")
            return

    def irSpectrum(self):
        """Generate a IR spectrum plot in a png file.

        It generates <jobname>_IR.png file.

        """
        frequencies = np.abs(self.freq[0])**0.5 / (2*self._pi*self._c)
        normodes = self.freq[1]
        intensity = self.intensity
        
        x = np.arange(500, 4000, 10)
        y = np.zeros(len(x))
        w = 5 #HWHM parameter
        for mode in range(len(frequencies)):
            y = y + ((1/self._pi) * intensity[mode] * w / ((x - frequencies[mode])**2 + w**2))

        y = 1.2*np.amax(y) - y
        fig, ax = plt.subplots()
        ax.set_xlabel('Frequencies [cm**-1]')
        ax.set_ylabel('Intensity')
        ax.plot(x,y)
        plt.savefig(f"{self.job}_IR.png")
        plt.show()
    

def xyzInp(job: str):
    """Import <jobname>.xyz file.

    It reads <jobname>.xyz file and store its information
    in coord and atmlist.

    Args:
        job (String): jobname, the name of molecule you calculate.

    Returns:
        list: atmlist, A list of atoms (String).
        numpy.ndarray: coord, An array of atom coordinates [Angstroms].

    """
    try:
        atmlist = []
        inpCoord = []
        f = open(job+'.xyz', 'r')
        datalist = f.readlines()
        for data in datalist:
            row = data.split()
        
            if len(row) != 4 :
                continue
        
            atm = row[0]
            xyz = [row[1], row[2], row[3]]
            atmlist.append(str(atm))
            inpCoord.extend(xyz)
        
        coord = np.array(inpCoord, dtype=float)
        return atmlist, coord
    
    except FileNotFoundError:
        print("There's no " + job +".xyz file. Please check the file name.")
        exit()


def main():
    job = input("Enter the job name: ")
    inp = xyzInp(job)
    atmlist = inp[0]
    coord = inp[1]
    print("Hessian Calculation Processing ...")
    hess = Hessian(coord, atmlist).hessian()
    harmo = Harmonic(hess, atmlist)
    eig = harmo.eig
    print(harmo)
    print("Dipole Moment Calculation Processing ...")
    ir = IR(eig, coord, atmlist)
    intensity = ir.xyzIntensity()
#    intensity = ir.normodeIntensity()
    print(ir)
    Plot(eig, job, intensity, atmlist)


if __name__ == "__main__":
    main()












