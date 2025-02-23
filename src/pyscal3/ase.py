import inspect


from ase.atoms import Atoms

import pyscal3.core as pc
from pyscal3.operations.centrosymmetry import calculate_centrosymmetry as _calculate_centrosymmetry
from pyscal3.operations.cna import calculate_cna as _calculate_cna


def _get_structure(ase_atoms: Atoms) -> pc.System:
    return pc.System(ase_atoms, format='ase')


def _get_updated_signature(signature: inspect.Signature) -> inspect.Signature:
    parameters = signature.parameters.copy()
    del parameters["system"]
    ase_atoms_parameter = inspect.Parameter(
        name="ase_atoms",
        kind=inspect.Parameter.POSITIONAL_OR_KEYWORD,
        annotation=Atoms,
    )
    parameters["ase_atoms"] = ase_atoms_parameter
    parameters.move_to_end("ase_atoms", last = False)
    return signature.replace(parameters=parameters.values())


def _wrap_function(input_function: callable) -> callable:
    sig_updated = _get_updated_signature(signature=inspect.signature(input_function))

    def ase_compatibility_wrapper(*args, **kwargs):
        input_dict = sig_updated.bind(*args, **kwargs).arguments
        input_dict["system"] = _get_structure(input_dict.pop("ase_atoms"))
        return input_function(**input_dict)

    ase_compatibility_wrapper.__name__ = input_function.__name__
    ase_compatibility_wrapper.__signature__ = sig_updated
    ase_compatibility_wrapper.__doc__ = input_function.__doc__.replace(
        "system: System object", "ase_atoms: ase.atoms.Atoms object"
    )
    return ase_compatibility_wrapper


calculate_centrosymmetry = _wrap_function(_calculate_centrosymmetry)
calculate_cna = _wrap_function(_calculate_cna)