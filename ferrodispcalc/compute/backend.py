from abc import ABC, abstractmethod

class ComputeBackend(ABC):
    @abstractmethod
    def get_averaged_structure(self, traj, select):
        pass 

    @abstractmethod
    def get_displacement(self, select, nl):
        pass

    @abstractmethod
    def get_polarization():
        pass