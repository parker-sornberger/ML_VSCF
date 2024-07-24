from ParseFequencyOutput import ParseFreq
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple

class ParseOutputVSCF(ParseFreq):
    def __init__(self, file, uses_scheme_1=False):
        super().__init__(file, supercell=None,upper_coord_block=False,sig_figs = None)
        self.scheme_1 = uses_scheme_1
        self._modes = set()
        self._md = {}
        self._step = 0
        self.bohr_steps = {}
        self._normal_step_factory = namedtuple("ClassicalAmplitudeSteps", 
                                                     ["Starting", "Ending", "Step"])
        self.classical_amplitude_steps = None
        
        self.reduced_masses = {}
        self._n3i = {}
        self._n4i = {}
        self._n2i1j = {}
        self._n3i1j = {}
        self._n2i2j = {}
        self.num_configs = 0
        self.convergence_energies = {}
        self.transition_energies = {}
        self._scheme_1_PES_arrays = {}
        self._zero = 0
        self._parse_fcs()
        
        
        
        num_modes = len(self._modes)

        one_d_shape = (num_modes, )
        two_d_shape = (num_modes, num_modes)
        self.n_2i = np.zeros(one_d_shape)
        self.n_3i = np.zeros(one_d_shape)
        self.n_4i = np.zeros(one_d_shape)
        self.n_2i_1j = np.zeros(two_d_shape)
        self.n_3i_1j = np.zeros(two_d_shape)
        self.n_2i_2j = np.zeros(two_d_shape)
        self._modes = sorted(self._modes)
        self.used_modes = self._modes
        self._construct_force_mats()
        self._calc_2nd_order()
        self.num_configs = len(self.convergence_energies)
    def _calc_2nd_order(self):
        wavs = np.array([self.modes[x]["WaveNumber"] for x in self.modes if x in self.used_modes])
        self.n_2i = (((wavs/(27.211386024367243*8065.73))**2)* 1822.8848)*0.5
    def _get_reduced_mass(self, line):
        mode, reduced_mass = line.strip().split()[-2:]
        self.reduced_masses[int(mode)] = float(reduced_mass)
    def _get_amp_step(self,line):
        line = line.strip().split()
        start = int(line[2])
        end = int(line[5])
        step = float(line[-1])
        self.classical_amplitude_steps = self._normal_step_factory(start, end, step)
    def _paired_bohr_step(self, line, modeline):
        
        modes = modeline.strip().split(":")[1].split()[:2]
        modes = tuple(map(int, modes))
        steps = line.strip().split()[1:3]
        steps = tuple(map(float, steps))
        self.bohr_steps[modes] = steps
        return True
        # if modes not in self.bohr_steps:
        #     self.bohr_steps[modes] = [steps]
        # else:
        #     self.bohr_steps[modes].append(steps)
    def _get_bohr_step(self, line, mode_index):
        modeline = self.lines[mode_index]
        if len(line.split())> 3:
            return self._paired_bohr_step(line, modeline)
        mode = int(modeline.strip().split(":")[1].split()[0])
        step = float(line.strip().split()[1])
        self.bohr_steps[mode] = step
        return False
        
        
    def _get_mode_tuple(self, line):
        line = line.replace(",", "").split(":")[-1].split("=")[0].replace("ETA(", "").replace(")", "").split()
        modetup = tuple(map(int, line))
        ms = set(modetup)
        for mode in ms:
            self._modes.add(mode)
        return modetup
    def _get_eta(self, line, modetup):
        eta = float(line.strip().split()[-2])
        if len(set(modetup)) == 1:
            if len(modetup) < 4:
                self._n3i[modetup] = eta
                
            else:
                self._n4i[modetup] = eta
            return
        if len(modetup) == 3:
            self._n2i1j[modetup] = eta
            return
        
        mode1, mode2 = sorted(set(modetup))
        m1_c = modetup.count(mode1)
        m2_c = modetup.count(mode2)
        if m1_c == m2_c:
            self._n2i2j[modetup] = eta
        else:
            self._n3i1j[modetup] = eta
        
        
        
    def _handle_eta_line(self, line):
        modetup = self._get_mode_tuple(line)
        eta = self._get_eta(line, modetup)
    def _get_config(self, line):
        config = int(line.strip().split(":")[-1])
        return config
    def _get_conv_energy(self, line, config):
        energy = float(line.strip().split()[-2])
        self.convergence_energies[config] = energy
    def _get_transition_energy(self, line, config):
        energy = float(line.strip().split()[-2])
        self.transition_energies[config] = energy
    
    def _get_paired_energies(self, chunk, modeline):
        modes = modeline.strip().split(":")[1].split()[:2]
        modes = tuple(map(int, modes))
        steps_x = []
        steps_y = []
        total = []
        harmonic = []
        delta = []
        #C    0.450   0.450  -9.167783E+02  -9.167797E+02     12  0.4E-11
        for line in chunk:


            _, stepx, stepy, tot, harm, _, delt = line.strip().split()
            delta.append(float(delt))

            
            steps_x.append(float(stepx))
            steps_y.append(float(stepy))
            total.append(float(tot))
            harmonic.append(float(harm))
        # harmonic = [x + self._zero for x in harmonic]
        # total = [x + self._zero for x in total]
        self._scheme_1_PES_arrays[modes] = dict(Range=(steps_x, steps_y), 
                                                Harmonic=harmonic, Total=total, DeltaE=delta)
    def _get_scheme1_PES(self, derivative_index, last_step_index, paired=False):
        
        step_ahead = 5 if not paired else 2
        steps_behind_to_modes = 2
        chunk = [line.strip() for line in self.lines[last_step_index + step_ahead:derivative_index] if line.strip()]
        modeline = self.lines[last_step_index - steps_behind_to_modes]
        if paired:
            return self._get_paired_energies(chunk, modeline)
        
        mode = int(modeline.strip().split(":")[1].split()[0])
        # C   -0.900   0.565213E-02   0.556310E-02     12 -0.2E-12
        #      0.000  -0.916781E+03  -0.916781E+03                         CENTRAL POINT
        steps = []
        total = []
        harmonic = []
        delta = []
        for line in chunk:

            if "CENTRAL" not in line:
                _, step, tot, harm, _, delt = line.strip().split()
                delta.append(float(delt))
            else:
                step, tot, harm, *_, = line.strip().split()
                self._zero = float(tot)
                tot = 0
                harm = 0
                delta.append(0)
            
            steps.append(float(step))
            total.append(float(tot))
            harmonic.append(float(harm))
            #"Harmonic", "Total", "DeltaE"
        harmonic = [x + self._zero for x in harmonic]
        total = [x + self._zero for x in total]
        self._scheme_1_PES_arrays[mode] = dict(Range=steps, Harmonic=harmonic, Total=total, DeltaE=delta)
        
        
        
    def _parse_fcs(self):
        reduced_mass_key = "REDUCED MASS OF MODE (AMU)"
        normalized_amp_step_keys = ("STARTING POINT:","ENDING POINT:", "STEP:")
        bohr_step_keys = ("STEP",  "BOHR") #covers both cases where 1 mode or 2 mode, 
        general_eta_key = "ETA("
        misspelled_config_key = "CONFIGUTARION NUMBER" #a misspelling that appears in some CRYSTAL versions
        config_key = "CONFIGURATION NUMBER"
        convergence_key = "CONVERGENCE ACHIEVED AT CYCLE"
        transition_key = "FUNDAMENTAL TRANSITION ENERGY" 
        
        potential_derivative_key = "POTENTIAL DERIVATIVES" #only matters to plot scheme 1 PES
        last_step_index = 0
        moved_to_paired_steps = False
        
        last_config = -1 
    
        for i, line in enumerate(self.lines):
            if reduced_mass_key in line:
                self._get_reduced_mass(line)
            if all(key in line for key in normalized_amp_step_keys):
                self._get_amp_step(line)
            if all(key in line for key in bohr_step_keys):
                moved_to_paired_steps = self._get_bohr_step(line, i-2)
                last_step_index = i
                
            if general_eta_key in line:
                self._handle_eta_line(line)
            if config_key in line or misspelled_config_key in line:
                last_config = self._get_config(line)
            if convergence_key in line:
                self._get_conv_energy(line, last_config)
            if transition_key in line:
                self._get_transition_energy(line, last_config)
            # if not self.scheme_1:
            #     continue
            # if potential_derivative_key in line:
            #     self._get_scheme1_PES(i, last_step_index, paired = moved_to_paired_steps)
                
    def _construct_force_mats(self):
        for key, eta in self._n3i.items():
            index = self._modes.index(key[0])
            self.n_3i[index] = eta
            
        for key, eta in self._n4i.items():
            index = self._modes.index(key[0])
            self.n_4i[index] = eta
        for key, eta in self._n2i1j.items():
            mode1, mode2 = sorted(set(key))
            m1_ind = self._modes.index(mode1)
            m2_ind = self._modes.index(mode2)
            self.n_2i_1j[(m1_ind, m2_ind)] = eta
        for key, eta in self._n3i1j.items():
            mode1, mode2 = sorted(set(key))
            m1_ind = self._modes.index(mode1)
            m2_ind = self._modes.index(mode2)
            self.n_3i_1j[(m1_ind, m2_ind)] = eta
        for key, eta in self._n2i2j.items():
            mode1, mode2 = sorted(set(key))
            m1_ind = self._modes.index(mode1)
            m2_ind = self._modes.index(mode2)
            self.n_2i_2j[(m1_ind, m2_ind)] = eta
        

    def _plot_2D(self, mode, which, zero, savefig_kwargs ):
        x = self._scheme_1_PES_arrays[mode]["Range"]
        y = self._scheme_1_PES_arrays[mode][which]
        plt.xlabel("Displacement Scale")
        plt.ylabel("Energy (A.U.)")
        if zero:
            y = [x - self._zero for x in y]

        plt.title(F"Mode {mode} {len(y)}-point PES")
        
        plt.plot(x, y)
        
        if isinstance(savefig_kwargs, dict):
            plt.savefig(**savefig_kwargs)
        plt.show()
    def _plot_3D(self, modes, which, zero, savefig_kwargs):


        xr, yr = self._scheme_1_PES_arrays[modes]["Range"]
        
        y = self._scheme_1_PES_arrays[modes][which]
        if zero:
            
            y = [x - self._zero for x in y]
            print(min(y))
        # _sortinds = np.argsort(xr)
        # xr = [xr[i] for i in _sortinds]
        # yr = [yr[i] for i in _sortinds]
        # y = [y[i] for i in _sortinds]
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(6,6))
        ax.azim = 60
        surf = ax.plot_trisurf(xr, yr, y, cmap='Blues')
        if isinstance(savefig_kwargs, dict):
            plt.savefig(**savefig_kwargs)
        plt.show()
    def plot_PES(self, mode, which = "Harmonic", zero=False, savefig_kwargs=None):
        
        if not self.scheme_1:
            print("Can only plot confirmations from Scheme 1 to calculate force constants...")
            return
        if which not in {"Harmonic", "Total", "DeltaE"}:
            return #raise sometrhint proban;y 
        if isinstance(mode, (tuple, list)):
            return self._plot_3D(mode, which, zero,savefig_kwargs)
        if mode not in self._scheme_1_PES_arrays:
            print(F"Don't have mode {mode} but we have {', '.join(str(x) for x in self._scheme_1_PES_arrays)}")
            return 
        self._plot_2D(mode, which, zero,savefig_kwargs)
        
    def write_VIBPOT(self,output_dir="", second_order_too=True, uncoupled=False):
        path = os.path.join(output_dir, "VIBPOT.DAT") if output_dir else "VIBPOT.DAT"
        with open(path, "w") as file:
            num_modes = len(self.used_modes)
            for i in range(num_modes):
                file.write(F"{float(self.eta_2i[i])}\n") if second_order_too else None
                file.write(F"{float(self.eta_3i[i])}\n")
                file.write(F"{float(self.eta_4i[i])}\n")
            if uncoupled:
                return
            for i in range(num_modes):
                for j in range(i+1, num_modes):
                    file.write(F"{float(self.eta_2i1j[(i,j)])}\n")
                    file.write(F"{float(self.eta_2i1j[(j,i)])}\n")
                    file.write(F"{float(self.eta_3i1j[(i,j)])}\n")
                    file.write(F"{float(self.eta_3i1j[(j,i)])}\n")
                    file.write(F"{float(self.eta_2i2j[(i,j)])}\n")
                    
            
    @property
    def eta_2i(self):
        return self.n_2i
    @property
    def eta_3i(self):
        return self.n_3i
    @property
    def eta_4i(self):
        return self.n_4i
    @property
    def eta_2i1j(self):
        return self.n_2i_1j
    @property
    def eta_3i1j(self):
        return self.n_3i_1j
    @property
    def eta_2i2j(self):
        return self.n_2i_2j
    @property
    def force_constants(self):

        return self.n_2i, self.n_3i, self.n_4i, self.n_2i_1j, self.n_3i_1j, self.n_2i_2j
       
if __name__=="__main__":
    import os
    path = r"my_totally_real_path"
    #file = "output_fully_coupled_VCI43_B3LYPD.out"
    #file = "my_cool_VSCF_calculation_file.out"
    file = "uncoupled_4_through_36_vscf.out"
    outfile = os.path.join(path, file)
    
    v = ParseOutputVSCF(outfile, uses_scheme_1=False)

















