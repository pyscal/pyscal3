"""
Create an annotated representation of the the system using pyscal_rdf
"""

def serialize(system, return_type='dict', outputfile=None):
	"""
	Serialize a system into a validated and annotated json representation.
	Needs pyscal-rdf module installed.

	Parameters
	----------
	return_type: Optional, {'dict', 'json', 'model', 'file'}
		if dict: a dictionary is returned,
		if json: a json dump is returned,
		if model: a pydantic class is returned,
		if file: output is written as a json file, outputfile needs
		to be specified

	outputfile: filename in which the system is to be serialised to

	Returns
	-------
	schema: output, type depends on `return_type`, see above.

	Note
	----
	Needs `pyscal_rdf`>=0.0.19 installation.
	Atom positions, types, and other properties are not serialised
	currently.
	"""
	try:
		from pyscal_rdf.schema import Sample
		import pyscal_rdf.properties as prp
	except ImportError:
		print("serialization needs pyscal_rdf module")

	#create sample
	sample = Sample()
	sample.material.element_ratio = system.composition

	if system._structure_dict is not None:
		sample.material.crystal_structure.altname = system.atoms._lattice
		spg = prp.get_space_group(system)
		sample.material.crystal_structure.spacegroup_symbol = spg[0]
		sample.material.crystal_structure.spacegroup_number = spg[1]
		sample.material.crystal_structure.unit_cell.bravais_lattice = prp.get_bravais_lattice(system)
		sample.material.crystal_structure.unit_cell.lattice_parameter = system.atoms._lattice_constant
		sample.material.crystal_structure.unit_cell.angle = list([prp.get_angle(system.box[0], system.box[1]),
															prp.get_angle(system.box[1], system.box[2]),
															prp.get_angle(system.box[2], system.box[0])])

	sample.simulation_cell.volume = system.volume
	sample.simulation_cell.number_of_atoms = system.natoms
	sample.simulation_cell.length = list(system.box_dimensions)
	sample.simulation_cell.vector = [list(x) for x in system.box]
	sample.simulation_cell.angle = list([prp.get_angle(system.box[0], system.box[1]),
															prp.get_angle(system.box[1], system.box[2]),
															prp.get_angle(system.box[2], system.box[0])])

	if return_type=='dict':
		return sample.dict()
	elif return_type=='json':
		return sample.model_dump_json()
	elif return_type=='model':
		return sample
	else:
		if outputfile is not None:
			with open(outputfile, 'w') as fout:
			    fout.write(sample.model_dump_json())


