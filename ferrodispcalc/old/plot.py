import matplotlib.pyplot as plt 
import numpy as np
from matplotlib import cm 
from glob import glob 

def calAngle(dx, dy):
    pp = np.sqrt(dx * dx + dy * dy)
    angle = np.arccos(dx / pp) / np.pi * 180.0
    index = np.where(dy < 0.0)
    angle[index] = 360.0 - angle[index]
    return angle

class Polter:
    def __init__(self,
                 disp_file: str,
                 origin_file: str) -> None:
        self.disp: np.ndarray = np.load(disp_file)
        self.origin: np.ndarray = np.load(origin_file)
        self.nframes = self.disp.shape[0]
        self.natoms = self.disp.shape[1]
        if self.nframes != self.origin.shape[0] or self.natoms != self.origin.shape[1]:
            raise ValueError("The shape of displacement and origin files do not match.")
        print("========Info========")
        print("Displacement file: ", disp_file)
        print("Origin file: ", origin_file)
        print("Number of atoms: ", self.natoms)
        print("Number of frames: ", self.nframes)
        print("====================")
    
    def _get_layer(self, origin: np.ndarray, disp: np.ndarray, key: int, tolerance: float=1) -> dict:
        layer_dict = {}
        round_d = np.floor(origin[:,key])

        for i in range(origin.shape[0]):
            layer_key = round_d[i]

            # Find the closest existing layer within the tolerance
            closest_layer = min(layer_dict.keys(), key=lambda x: abs(x-layer_key) if abs(x-layer_key) <= tolerance else float('inf'), default=layer_key)

            if abs(closest_layer - layer_key) <= tolerance:
                layer_key = closest_layer  # Consolidate to the closest layer

            if layer_key not in layer_dict:
                layer_dict[layer_key] = []
            
            layer_dict[layer_key].append([origin[i], disp[i]])

        return layer_dict

    def plot_average_layer_polarization_mxs(self, timestep: int=250,size: list[int]=[10,10,10],
                                            save: bool=False, prefix: str=None) -> None:
        disp = self.disp.copy()[timestep:]
        disp = np.mean(disp, axis=0)
        dx, dy, dz = disp[:,0], disp[:,1], disp[:,2]
        dx = dx.reshape(size)
        dy = dy.reshape(size)
        dz = dz.reshape(size)
        for i in range(3):
            for j in range(size[i]):
                if i==0: #yz
                    dy0 = dy[j,:,:]
                    dz0 = dz[j,:,:]
                    angle = calAngle(dy0.flatten(), dz0.flatten())
                    profile = angle.reshape(dy0.shape)
                    plt.quiver(dy0, dz0)
                    plt.imshow(profile, cmap=cm.hsv, vmax=360, vmin=0, aspect=1.0, origin='lower')
                    plt.xlabel("Y")
                    plt.ylabel("Z")
                    plt.title("YZ-Plane")
                    plt.savefig(f'{prefix}-YZ-Plane-Layer{j}.png',bbox_inches='tight',dpi=300)
                    plt.close()
                elif i==1: #xz
                    dx0 = dx[:,j,:]
                    dz0 = dz[:,j,:]
                    angle = calAngle(dx0.flatten(), dz0.flatten())
                    profile = angle.reshape(dx0.shape)
                    plt.quiver(dx0, dz0)
                    plt.imshow(profile, cmap=cm.hsv, vmax=360, vmin=0, aspect=1.0, origin='lower')
                    plt.xlabel("X")
                    plt.ylabel("Z")
                    plt.title("XZ-Plane")
                    plt.savefig(f'{prefix}-XZ-Plane-Layer{j}.png',bbox_inches='tight',dpi=300)
                    plt.close()
                elif i==2: #xy
                    dx0 = dx[:,:,j]
                    dy0 = dy[:,:,j]
                    angle = calAngle(dx0.flatten(), dy0.flatten())
                    profile = angle.reshape(dx0.shape)
                    plt.quiver(dx0, dy0)
                    plt.imshow(profile, cmap=cm.hsv, vmax=360, vmin=0, aspect=1.0, origin='lower')
                    plt.xlabel("X")
                    plt.ylabel("Y")
                    plt.title("XY-Plane")
                    plt.savefig(f'{prefix}-XY-Plane-Layer{j}.png',bbox_inches='tight',dpi=300)
                    plt.close()
            

        
    def plot_average_layer_polarization(self, timestep: int=250,save: bool=False, prefix: str=None) -> None:

        # get the average displacement and origin
        disp = self.disp.copy()[timestep:]
        origin = self.origin.copy()[timestep:]
        disp = np.mean(disp, axis=0)
        origin = np.mean(origin, axis=0)
        
        for i in range(3):
            layer_dict = self._get_layer(origin, disp, i)
            for key in layer_dict:
                origin_list, disp_list = zip(*layer_dict[key])
                origin_list = np.array(origin_list)
                disp_list = np.array(disp_list)

                if i==0: # yz plane
                    x,y = origin_list[:,1], origin_list[:,2]
                    dx,dy = disp_list[:,1], disp_list[:,2]
                    tit = 'YZ-Plane'
                    xlabel="Y"
                    ylabel="Z"
                elif i==1: # xz plane
                    x,y = origin_list[:,0], origin_list[:,2]
                    dx,dy = disp_list[:,0], disp_list[:,2]
                    tit = 'XZ-Plane'
                    xlabel="X"
                    ylabel="Z"
                elif i==2: # xy plane
                    x,y = origin_list[:,0], origin_list[:,1]
                    dx,dy = disp_list[:,0], disp_list[:,1]
                    tit = 'XY-Plane'
                    xlabel="X"
                    ylabel="Y"
                
                angles = np.arctan2(dy,dx)
                colormap = cm.hsv
                colors = colormap(angles)
                plt.quiver(x,y,dx,dy,color=colors)
                plt.title(tit)
                plt.xlabel(xlabel)
                plt.ylabel(ylabel)
                if save:
                    plt.savefig(f'{prefix}-{tit}-{key}.png',dpi=600)
                plt.close()
    def plot_displacement_distribution(self, prefix: str, save: bool=True,
                                       tit: str=None):
        disp = self.disp.reshape(-1,3)
        d = np.linalg.norm(disp,axis=1)
        
        # plot dx
        hist,bin_edges = np.histogram(disp[:,0], bins=100, density=True)
        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        plt.plot(bin_centers, hist, linestyle='-', linewidth=2, label=f'dx')

        # plot dy
        hist,bin_edges = np.histogram(disp[:,1], bins=100, density=True)
        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        plt.plot(bin_centers, hist, linestyle='-', linewidth=2, label=f'dy')

        # plot dz
        hist,bin_edges = np.histogram(disp[:,2], bins=100, density=True)
        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        plt.plot(bin_centers, hist, linestyle='-', linewidth=2, label=f'dz')

        # plot d
        hist,bin_edges = np.histogram(d, bins=100, density=True)
        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        plt.plot(bin_centers, hist, linestyle='-', linewidth=2, label=f'd')

        # config
        plt.xlabel('Displacement (Å)')
        plt.ylabel('Frequency')
        plt.legend()
        plt.grid()
        if save:
            plt.savefig(f"{prefix}disp_distribution.png",dpi=600)
        plt.close()
    
    def plot_phase_percentage(self, file_list: list[str], save: bool=True, 
                              tit: str=None, prefix: str=None):
        phase_percentages = {'C': [], 'R': [], 'O': [], 'T': []}
        R,O,T,C=0,0,0,0

        for file in file_list:
            d: np.ndarray = np.load(file)
            R,O,T,C=0,0,0,0

            # calculate displacement of each atom
            for frame in d:
                for atom in frame:
                    d = np.linalg.norm(atom)
                    if d == 0:
                        continue
                    if d<0.1:
                        C+=1
                    else:
                        polar_count = sum(abs(k) > d/np.sqrt(6) for k in atom)
                        if polar_count == 3:
                            R+=1
                        elif polar_count == 1:
                            T+=1
                        elif polar_count == 2:
                            O+=1
            total = R+O+T+C
            phase_percentages['C'].append(C / total * 100)
            phase_percentages['R'].append(R / total * 100)
            phase_percentages['O'].append(O / total * 100)
            phase_percentages['T'].append(T / total * 100)
        bar_width = 0.2
        r1 = np.arange(len(file_list))
        r2 = [x + bar_width for x in r1]
        r3 = [x + bar_width for x in r2]
        r4 = [x + bar_width for x in r3]
        plt.bar(r2, phase_percentages['R'], color='r', width=bar_width, label='R')
        plt.bar(r3, phase_percentages['O'], color='g', width=bar_width, label='O')
        plt.bar(r4, phase_percentages['T'], color='y', width=bar_width, label='T')
        plt.bar(r1, phase_percentages['C'], color='b', width=bar_width, label='C')
        print(phase_percentages)
        plt.xlabel('Temperature', fontweight='bold')
        plt.xticks([r + bar_width for r in range(len(file_list))], file_list)
        plt.ylabel('Phase Percentage', fontweight='bold')
        plt.title(tit)
        plt.legend()
        if save:
            plt.savefig('phase_percentage.png', dpi=600)
        plt.close()







a = Polter('0.3750-origin.npy','0.3750-disp.npy')
file_list = glob('*-disp.npy')
print(file_list)
a.plot_phase_percentage(file_list=file_list,save=True,
                        tit='Precentage', prefix='t')