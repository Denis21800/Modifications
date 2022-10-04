import os.path
from pathlib import Path

import pyrosetta

from logger.tasklogger import TaskLogger
from pdb_tools.clean_pdb import CleanPDB
from config.config import Config
from Bio.PDB.Polypeptide import aa3

from modifiations_tools.const import NC_AAS_PATCH
from modifiations_tools.mutation_task import MutationTask
from transformers.optimization import Minimization
from modifiations_tools.scores import ScoreMetrics
from visualization.visualize import VisualizeMng


def rosetta_configure():
    pyrosetta.init(extra_options="-mute core -mute basic")


class MutationMngr(object):
    def __init__(self):
        self.config = None
        self.pdb_file = None
        self.residues_info = None
        self.logger = None
        self.task_id = -1
        self.mutation_map_info = {}

    # mutate_residues: [(positions, new res)]
    def create_task(self,
                    pdb_path,
                    model_id,
                    chain_id,
                    mutate_residues
                    ):
        self.task_id = -1
        try:
            self.config = Config()
            assert os.path.exists(pdb_path), f"Can't find input pdb: {pdb_path}"
            self.task_id = self.config.make_working_space_dir()
            logger_params = self.config.get_logging_params()
            self.logger = TaskLogger(log_params=logger_params)
            tmp_dir = self.config.get_tmp_dir()
            clean_pdb = CleanPDB(input_pdb=pdb_path,
                                 model_id=model_id,
                                 chain_id=chain_id,
                                 out_pdb_path=tmp_dir)
            self.pdb_file, self.residues_info = clean_pdb.clean()
            assert os.path.exists(self.pdb_file), f"Can't read input pdb: {pdb_path}"
            assert self.residues_info, f"Invalid pdb structure: {pdb_path}"
            check_result = self.__check_inputs(mutate_residues=mutate_residues)
            assert check_result, f"Invalid pdb structure: {pdb_path}"
        except Exception as e:
            return -1, str(e)
        return self.task_id, None

    def mutate(self, task_id):
        results = False
        try:
            assert task_id == self.task_id and task_id != -1, f"Invalid task id {task_id}"
            self.__log_event(self.logger.INFO, f"Processing task: {task_id}")
            mutation_task = MutationTask(reference_pdb_path=self.pdb_file,
                                         mutation_info=self.mutation_map_info,
                                         config_params=self.config.get_mutation_params()
                                         )
            self.__log_event(self.logger.INFO, str(mutation_task.reference_pose))
            repack_result = mutation_task.mutate()
            for key in repack_result:
                self.__log_event(self.logger.INFO, f"Repacking: {key} \n {repack_result[key]}")
            self.__log_event(self.logger.INFO, "Relaxation started...")
            min_packer = Minimization(score_fn=mutation_task.score_fn)
            transformed_reference_pose = min_packer.minimize_std(mutation_task.reference_pose)
            transformed_target_pose = min_packer.minimize_std(mutation_task.target_pose)
            poses_data = [mutation_task.reference_pose, mutation_task.target_pose,
                          transformed_reference_pose, transformed_target_pose]

            input_pdb_id = Path(self.pdb_file).stem
            pdb_out_params = self.config.get_pdb_out_params()
            pdb_out_dir = pdb_out_params['pdb_out_dir']
            assert os.path.exists(pdb_out_dir), f"Can't access to pdb out folder {pdb_out_dir}"
            for i, pose in enumerate(poses_data):
                postfix = pdb_out_params[list(pdb_out_params.keys())[i + 1]]
                out_file = pdb_out_dir / f"{input_pdb_id}_{postfix}.pdb"
                pyrosetta.dump_pdb(pose, str(out_file))
            scores = []
            scores_data_path = self.config.get_results_data_params()

            for i, pose in enumerate(poses_data):
                score = ScoreMetrics(
                    pose=pose,
                    score_fn=mutation_task.score_fn,
                    mutation_points=self.mutation_map_info
                )
                score.calculate()
                data_file = scores_data_path[list(scores_data_path.keys())[i]]
                score.dump(data_file)
                scores.append(score)
            try:
                vis_params = self.config.get_visualisation_params()
                if vis_params['save']:
                    score_results = [sc.results for sc in scores]
                    vis_mng = VisualizeMng(vis_params)
                    vis_mng.init_data(
                        reference=score_results[0],
                        target=score_results[1],
                        reference_transformed=score_results[2],
                        target_transformed=score_results[3]
                    )
                    vis_mng.plot()
            except Exception as e:
                self.__log_event(self.logger.ERROR, f'Visualization error: {str(e)}')
            self.__log_event(self.logger.INFO, f"Task {self.task_id} has been completed successfully")
            results = True
        except Exception as e:
            if not self.logger:
                raise
            self.__log_event(self.logger.CRITICAL, f"Undefined exception': {str(e)}")
        finally:
            self.config.clear_tmp_ws()

        return results

    def __check_inputs(self, mutate_residues):
        if not self.residues_info:
            self.__log_event(self.logger.ERROR, "Error: Invalid input file: no residues found")

        for item in mutate_residues:
            position, new_res_code = item
            mutate_record = self.residues_info.get(position)
            if not mutate_record:
                self.__log_event(self.logger.WARNING, f"Error: can't find residues {position}")
            reference_res = mutate_record['res_aa']
            order_pos = mutate_record['res_order']
            if new_res_code not in list(aa3) and new_res_code not in NC_AAS_PATCH:
                self.__log_event(self.logger.WARNING, f"ignore unrecognized code: {new_res_code}")
            self.__log_event(self.logger.INFO, f"Modify: {reference_res}->{new_res_code}")
            self.mutation_map_info.update({order_pos: new_res_code})

        if not self.mutation_map_info:
            self.__log_event(self.logger.ERROR, "No valid positions found")
            return False
        return True

    def __log_event(self, level, message):
        if self.logger:
            self.logger.log_event(level, message)
