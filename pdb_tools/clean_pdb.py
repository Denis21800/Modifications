import os.path

from Bio.PDB import PDBParser, PDBIO, Select
from pyrosetta.toolbox import cleanATOM


class CleanPDB(object):
    def __init__(self,
                 input_pdb,
                 model_id,
                 chain_id,
                 out_pdb_path):
        assert os.path.exists(input_pdb)
        self.input_pdb = input_pdb
        self.out_pdb_path = out_pdb_path / os.path.basename(input_pdb)
        self.model_id = model_id
        self.chain_id = chain_id
        self.residues_map = {}

    def clean(self):
        cleanATOM(str(self.input_pdb), out_file=self.out_pdb_path)
        assert os.path.exists(self.out_pdb_path)
        self.extract_chain()
        return self.out_pdb_path, self.residues_map

    def extract_chain(self):
        parser = PDBParser(PERMISSIVE=1)
        structure = parser.get_structure(0, file=self.out_pdb_path)
        for model in structure:
            if model.get_id() != self.model_id:
                continue
            for chain in model:
                chain_id = chain.get_id()
                if chain_id != self.chain_id:
                    continue
                for idx, res in enumerate(chain.get_residues()):
                    _, res_id, _ = res.get_id()
                    res_name = res.get_resname()
                    self.residues_map.update({res_id: {'chain': chain_id,
                                                       'res_aa': res_name,
                                                       'res_order': idx + 1}
                                              })

        pdb_io = PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save(str(self.out_pdb_path),
                    ChainSelect(chain_id=self.chain_id,
                                model_id=self.model_id
                                ))


class ChainSelect(Select):
    def __init__(self,
                 model_id,
                 chain_id,
                 ):
        self.selected_model = model_id
        self.selected_chain = chain_id

    def accept_model(self, model):
        return model.get_id() == self.selected_model

    def accept_chain(self, chain):
        return chain.get_id() == self.selected_chain
