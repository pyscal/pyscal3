import inspect

from ase.atoms import Atoms

import pyscal3.core as pc
from pyscal3.operations.centrosymmetry import calculate_centrosymmetry as _calculate_centrosymmetry
from pyscal3.operations.calculations import calculate_q as _calculate_q
from pyscal3.operations.cna import (
    calculate_cna as _calculate_cna,
    identify_diamond as _identify_diamond,
)
from pyscal3.operations.identify import find_neighbors


def _get_structure(ase_atoms: Atoms) -> pc.System:
    return pc.System(ase_atoms, format='ase')


def _get_updated_signature(signature: inspect.Signature, add_find_neighbors: bool = False) -> inspect.Signature:
    parameters = signature.parameters.copy()
    del parameters["system"]
    parameters["ase_atoms"] = inspect.Parameter(
        name="ase_atoms",
        kind=inspect.Parameter.POSITIONAL_OR_KEYWORD,
        annotation=Atoms,
    )
    if add_find_neighbors:
        for k, v in inspect.signature(find_neighbors).parameters.items():
            if k != "system":
                parameters[k] = v
    parameters.move_to_end("ase_atoms", last = False)
    return signature.replace(parameters=parameters.values())


def _wrap_function(funct: callable, add_find_neighbors: bool = False) -> callable:
    sig_updated = _get_updated_signature(
        signature=inspect.signature(funct),
        add_find_neighbors=add_find_neighbors,
    )

    def ase_compatibility_wrapper(*args, **kwargs):
        input_dict = sig_updated.bind(*args, **kwargs).arguments
        pc_system = _get_structure(input_dict.pop("ase_atoms"))
        if add_find_neighbors:
            find_neighbor_keys = inspect.signature(find_neighbors).parameters.keys()
            find_neighbors_kwargs = {
                k: input_dict.pop(k)
                for k in find_neighbor_keys
                if k != "system" and k in input_dict
            }
            pc_system.find.neighbors(**find_neighbors_kwargs)
        input_dict["system"] = pc_system
        return funct(**input_dict)

    ase_compatibility_wrapper.__name__ = funct.__name__
    ase_compatibility_wrapper.__signature__ = sig_updated
    ase_compatibility_wrapper.__doc__ = funct.__doc__.replace(
        "system: System object", "ase_atoms: ase.atoms.Atoms object"
    )
    return ase_compatibility_wrapper


calculate_centrosymmetry = _wrap_function(funct=_calculate_centrosymmetry, add_find_neighbors=False)
calculate_cna = _wrap_function(funct=_calculate_cna, add_find_neighbors=False)
calculate_steinhardt_parameter = _wrap_function(funct=_calculate_q, add_find_neighbors=True)
calculate_diamond_structure = _wrap_function(funct=_identify_diamond, add_find_neighbors=False)