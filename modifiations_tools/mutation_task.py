import os.path
import sys

import pyrosetta
from pyrosetta import rosetta
from pyrosetta.rosetta.core import conformation
from pyrosetta.rosetta.core import chemical
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.toolbox.mutants import restrict_non_nbrs_from_repacking

from modifiations_tools.const import NC_AAS_PATCH


class MutationTask(object):
    def __init__(self,
                 reference_pdb_path,
                 mutation_info,
                 config_params,
                 ):
        assert os.path.exists(reference_pdb_path)
        self.reference_path = str(reference_pdb_path)
        self.mutation_info = mutation_info
        self.score_fn = pyrosetta.get_fa_scorefxn()
        self.reference_pose = None
        self.target_pose = None
        self.__init_structure()
        self.repack_radius = config_params['repack_radius']

    def __init_structure(self):
        self.reference_pose = pyrosetta.pose_from_pdb(self.reference_path)
        self.target_pose = pyrosetta.Pose()
        self.target_pose.assign(self.reference_pose)

    def mutate(self):
        repacks_log = {}
        for residue_position in self.mutation_info:
            self.__mutate_aa(residue_position, self.mutation_info[residue_position])
            repack_task = self.__repack_residue(residue_position)
            repacks_log.update({self.mutation_info[residue_position]: str(repack_task)})
        return repacks_log

    def __repack_residue(self, mutant_position):
        repack_task = pyrosetta.rosetta.core.pack.task.TaskFactory.create_packer_task(self.target_pose)
        repack_task = restrict_non_nbrs_from_repacking(self.target_pose,
                                                       mutant_position,
                                                       repack_task,
                                                       pack_radius=self.repack_radius)
        packer = PackRotamersMover(self.score_fn, repack_task)
        packer.apply(self.target_pose)
        return repack_task

    def __mutate_aa(self, residue_position, aa_name):
        residue = self.target_pose.residue(residue_position)

        if residue.name() in ('CYS:disulfide', 'CYD'):
            self.__di_sulfide_cys_repack(residue)

        chm_manager = chemical.ChemicalManager.get_instance()
        rts_factory = chm_manager.residue_type_set("fa_standard")
        if aa_name in NC_AAS_PATCH:
            cognate_res_type = rts_factory.name_map(NC_AAS_PATCH[aa_name]['cognateAA'])
            mutation_res_type = rts_factory.get_residue_type_with_variant_added(cognate_res_type,
                                                                                NC_AAS_PATCH[aa_name]['type'])
        else:
            mutation_res_type = rts_factory.name_map(aa_name)

        if residue_position == 1:
            try:
                variant_type = chemical.VariantType.LOWER_TERMINUS_VARIANT
                mutation_res_type = rts_factory.get_residue_type_with_variant_added(mutation_res_type,
                                                                                    variant_type)
            except AttributeError:
                variant_type = "LOWER_TERMINUS"
                mutation_res_type = rts_factory.get_residue_type_with_variant_added(mutation_res_type, variant_type)
        elif residue == self.target_pose.total_residue():
            try:
                variant_type = chemical.VariantType.UPPER_TERMINUS_VARIANT
                mutation_res_type = rts_factory.get_residue_type_with_variant_added(mutation_res_type, variant_type)
            except AttributeError:
                variant_type = "UPPER_TERMINUS"
                mutation_res_type = rts_factory.get_residue_type_with_variant_added(mutation_res_type, variant_type)
        mutate_residue = conformation.ResidueFactory.create_residue(mutation_res_type)
        self.target_pose.replace_residue(residue_position, mutate_residue, orient_backbone=True)

    def __di_sulfide_cys_repack(self, residue):
        try:
            di_sulfide_partner = residue.residue_connection_partner(
                residue.n_residue_connections())
        except AttributeError:
            di_sulfide_partner = residue.residue_connection_partner(
                residue.n_current_residue_connections())
        temp_pose = pyrosetta.Pose()
        temp_pose.assign(self.target_pose)
        conformation.change_cys_state(residue, 'CYS', temp_pose.conformation())
        conformation.change_cys_state(di_sulfide_partner, 'CYS', temp_pose.conformation())
        self.target_pose = temp_pose
