import matplotlib.pyplot as plt 
import numpy as np
from matplotlib import cm 
from glob import glob 
from matplotlib.animation import FuncAnimation
import os 

def calAngle(dx, dy):
    pp = np.sqrt(dx * dx + dy * dy)
    angle = np.arccos(dx / pp) / np.pi * 180.0
    index = np.where(dy < 0.0)
    angle[index] = 360.0 - angle[index]
    return angle


class DispPloter:
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

    def plot_average_layer_polarization_mxs(self, timestep: int=250,size: list[int]=[10,10,10],
                                            save: bool=False, prefix: str=None,
                                            dir: str=None) -> None:
        disp = self.disp.copy()[timestep:]
        disp = np.mean(disp, axis=0)
        dx, dy, dz = disp[:,0], disp[:,1], disp[:,2]
        dx = dx.reshape(size)
        dy = dy.reshape(size)
        dz = dz.reshape(size)

        if dir is not None:
            os.makedirs(dir, exist_ok=True)
        
        for i in range(3):
            for j in range(size[i]):
                if i==0: #yz
                    dy0 = dy[j,:,:]
                    dz0 = dz[j,:,:]
                    angle = calAngle(dy0.flatten(), dz0.flatten())
                    profile = angle.reshape(dy0.shape)
                    plt.quiver(dy0, dz0)
                    plt.imshow(profile, cmap=cm.hsv, vmax=360, vmin=0, aspect=1.0, origin='lower')
                    plt.title("YZ-Plane")
                    plt.savefig(f'{dir}/{prefix}-YZ-Plane-Layer{j}.png',bbox_inches='tight',dpi=300)
                    plt.close()
                elif i==1: #xz
                    dx0 = dx[:,j,:]
                    dz0 = dz[:,j,:]
                    angle = calAngle(dx0.flatten(), dz0.flatten())
                    profile = angle.reshape(dx0.shape)
                    plt.quiver(dx0, dz0)
                    plt.imshow(profile, cmap=cm.hsv, vmax=360, vmin=0, aspect=1.0, origin='lower')
                    plt.title("XZ-Plane")
                    plt.savefig(f'{dir}/{prefix}-XZ-Plane-Layer{j}.png',bbox_inches='tight',dpi=300)
                    plt.close()
                elif i==2: #xy
                    dx0 = dx[:,:,j]
                    dy0 = dy[:,:,j]
                    angle = calAngle(dx0.flatten(), dy0.flatten())
                    profile = angle.reshape(dx0.shape)
                    plt.quiver(dx0, dy0)
                    plt.imshow(profile, cmap=cm.hsv, vmax=360, vmin=0, aspect=1.0, origin='lower')
                    plt.title("XY-Plane")
                    plt.savefig(f'{dir}/{prefix}-XY-Plane-Layer{j}.png',bbox_inches='tight',dpi=300)
                    plt.close()
    def plot_polarization_time(self,plane: str, layer_index: int, size: list[int,int,int],save: bool=False, prefix: str=None, dir: str=None) -> None:
        fig, ax = plt.subplots()

        def update(frame):
            plt.cla()
            disp = self.disp.copy()[frame]
            dx, dy, dz = disp[:,0], disp[:,1], disp[:,2]
            dx = dx.reshape(size)
            dy = dy.reshape(size)
            dz = dz.reshape(size)
            if plane == 'yz':
                dy0 = dy[layer_index,:,:]
                dz0 = dz[layer_index,:,:]
                angle = calAngle(dy0.flatten(), dz0.flatten())
                profile = angle.reshape(dy0.shape)
                ax.quiver(dy0, dz0)
                ax.imshow(profile, cmap=cm.hsv, vmax=360, vmin=0, aspect=1.0, origin='lower')
                ax.set_title(f"YZ-Plane Layer {layer_index} - Frame {frame}")
            elif plane == 'xz':
                dx0 = dx[:,layer_index,:]
                dz0 = dz[:,layer_index,:]
                angle = calAngle(dx0.flatten(), dz0.flatten())
                profile = angle.reshape(dx0.shape)
                ax.quiver(dx0, dz0)
                ax.imshow(profile, cmap=cm.hsv, vmax=360, vmin=0, aspect=1.0, origin='lower')
                ax.set_title(f"XZ-Plane Layer {layer_index} - Frame {frame}")
            elif plane == 'xy':
                dx0 = dx[:,:,layer_index]
                dy0 = dy[:,:,layer_index]
                angle = calAngle(dx0.flatten(), dy0.flatten())
                profile = angle.reshape(dx0.shape)
                ax.quiver(dx0, dy0)
                ax.imshow(profile, cmap=cm.hsv, vmax=360, vmin=0, aspect=1.0, origin='lower')
                ax.set_title(f"XY-Plane Layer {layer_index} - Frame {frame}")

        ani = FuncAnimation(fig, update, frames=self.nframes, repeat=False,interval=50)
        plt.show()

        #if save:
        #    if dir is not None:
        #        os.makedirs(dir, exist_ok=True)
        #    ani.save(f'{dir}/{prefix}-{plane}-Layer{layer_index}.gif', writer='pillow', fps=10)
        #plt.close()
    def plot_layer_polarization(self, time: int=None, average: bool=True, start_idx: int=None, layer_tolerance: float=1.0):
        
        # get the original point
        origin = self.origin[0].copy()
        
        # get desired displacement
        if average:
            disp = np.mean(self.disp[start_idx:].copy(), axis=0)
        else:
            disp = self.disp[time].copy()
        
        # according to the origin's coord, find the layer
        x_layer_dict = {}  # key: x , value: list of atoms
        x_layer_key = []
        y_layer_dict = {}
        y_layer_key = []
        z_layer_dict = {}
        z_layer_key = []

        for idx, (pos, d) in enumerate(zip(origin, disp)):
            x, y, z = pos

            # yz
            found_layer= False
            for key in x_layer_dict:
                if abs(x-key) < layer_tolerance:
                    x_layer_dict[key].append((y,z,d[1], d[2]))
                    found_layer = True
                    break 
            
            if not found_layer:
                x_layer_key.append(x)
                x_layer_dict[x] = [(y,z,d[1], d[2])]
            
            # xz
            found_layer= False
            for key in y_layer_dict:
                if abs(y-key) < layer_tolerance:
                    y_layer_dict[key].append((x,z,d[0], d[2]))
                    found_layer = True
                    break 
            
            if not found_layer:
                y_layer_key.append(y)
                y_layer_dict[y] = [(x,z,d[0], d[2])]
            
            # xy
            found_layer= False
            for key in z_layer_dict:
                if abs(z-key) < layer_tolerance:
                    z_layer_dict[key].append((x,y,d[0], d[1]))
                    found_layer = True
                    break 
            
            if not found_layer:
                z_layer_key.append(z)
                z_layer_dict[z] = [(x,y,d[0], d[1])]
        print(z_layer_key)
        # plot xy 
        for i, key in enumerate(sorted(z_layer_key)):
            x,y,dx,dy = zip(*z_layer_dict[key])
            plt.quiver(x,y,dx,dy)
            plt.savefig(f'{key}.png')
            plt.close()

            