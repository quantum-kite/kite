import numpy as np

class Sublattice:
    def __init__(self, position, energy, alias_id):
        self.position = position
        self.energy = energy
        self.alias_id = alias_id

class HoppingFamily:
    def __init__(
        self,
        relative_index,
        from_sub: str,
        to_sub: str,
        energy,
        from_id: int,
        to_id: int,
    ):
        self.relative_index = relative_index
        self.from_sub = from_sub
        self.to_sub = to_sub
        self.energy = energy
        self.terms = [HoppingTerms(relative_index, from_id, to_id)]

class HoppingTerms:
    def __init__(self, relative_index: list, from_id: int, to_id: int):
        self.relative_index = relative_index
        self.from_id = from_id
        self.to_id = to_id

# Class that replaces the pybinding for lattice generation
class Lattice:
    def __init__(self, a1, a2=None, a3=None):
        self._vectors = [v for v in (a1, a2, a3) if v is not None]
        self.sublattices = {}
        self.hoppings = {}
        self._hopping_energies_map = {}
        self.nsub = 0
        self.nhop = 0

    @property
    def vectors(self) -> list:
        """Primitive lattice vectors"""
        return self._vectors

    def add_one_hopping(
        self, relative_index: list[int], from_sub: str, to_sub: str, hop_name_or_energy
    ) -> None:
        """Does not add the complex conjugate hopping"""
        name = f"Hop{self.nhop}"

        if type(hop_name_or_energy) == str:
            if hop_name_or_energy in self._hopping_energies_map:
                hop_val = self._hopping_energies_map[hop_name_or_energy]
            else:
                raise TypeError(f"The hopping label has not been defined")
        else:
            hop_val = complex(hop_name_or_energy)

        self.hoppings[name] = HoppingFamily(
            np.asarray(relative_index),
            from_sub,
            to_sub,
            np.array([[hop_val]]),
            self.sublattices[from_sub].alias_id,
            self.sublattices[to_sub].alias_id,
        )
        self.nhop += 1

    def add_hoppings(self, *hoppings) -> None:
        """Does not add the complex conjugate hopping"""
        for hop in hoppings:
            self.add_one_hopping(*hop)

    def add_one_sublattice(
        self, name: str, position: list[float], onsite_energy: float = 0.0
    ) -> None:
        if name in self.sublattices:
            raise ValueError(f"Sublattice '{name}' already exists")

        self.sublattices[name] = Sublattice(
            position, np.array([[onsite_energy]]), self.nsub
        )
        self.nsub += 1

    def add_sublattices(self, *sublattices) -> None:
        for sub in sublattices:
            self.add_one_sublattice(*sub)

    def register_hopping_energies(self, mapping: dict[str, float]) -> None:
        """Register a mapping of user-friendly names to hopping energies"""
        self._hopping_energies_map = mapping
