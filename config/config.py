import configparser
import os.path
import shutil
import uuid
from configparser import ConfigParser
from pathlib import Path


class Config(object):
    CONFIG_PATH = 'config.ini'

    def __init__(self):
        self.config = ConfigParser()
        assert os.path.exists(self.CONFIG_PATH), f"Can't find config file: {self.CONFIG_PATH}"
        self.config.read(self.CONFIG_PATH)
        self.working_dir = Path(self.config.get('working', 'working_dir'))
        tmp_dir = self.config.get('working', 'tmp_dir')
        assert os.path.exists(self.working_dir)
        self.space_id = uuid.uuid4()
        self.working_spase_dir = self.working_dir / str(self.space_id)
        self.tmp_dir = self.working_spase_dir / tmp_dir
        results_dir = self.config.get('working', 'out_dir')
        self.results_dir = self.working_spase_dir / results_dir
        vis_dir = self.config.get('visualization', 'out_dir')
        self.out_visualization_dir = self.results_dir / vis_dir
        self.json_data_dir = self.config.get('working', 'results_data_dir')
        self.pdb_out_dir = self.config.get('working', 'results_pdb_dir')
        self.json_data_dir = self.results_dir / self.json_data_dir
        self.pdb_out_dir = self.results_dir / self.pdb_out_dir
        self.log_dir = self.config.get('logging', 'log_dir')
        self.log_dir = self.working_spase_dir / self.log_dir

    def make_working_space_dir(self):
        self.working_spase_dir.mkdir(exist_ok=True)
        self.tmp_dir.mkdir(exist_ok=True)
        self.results_dir.mkdir(exist_ok=True)
        self.out_visualization_dir.mkdir(exist_ok=True)
        self.json_data_dir.mkdir(exist_ok=True)
        self.pdb_out_dir.mkdir(exist_ok=True)
        self.log_dir.mkdir(exist_ok=True)
        return self.space_id

    def get_vis_dir(self):
        return self.out_visualization_dir

    def get_tmp_dir(self):
        return self.tmp_dir

    def get_results_dir(self):
        return self.results_dir

    def get_logging_params(self):
        params = {}
        try:
            logging = self.config.get('logging', 'logging')
            params['logging'] = bool(logging)
            out = self.config.get('logging', 'out')
            params['out'] = out
            params['level'] = self.config.get('logging', 'level')
            if not logging or out == 'console':
                shutil.rmtree(self.log_dir)
                return params

            log_file = self.config.get('logging', 'log_file')
            log_file = self.log_dir / log_file
            params['log_file'] = str(log_file)

        except (configparser.NoSectionError, KeyError, ValueError):
            raise
            if os.path.exists(self.log_dir):
                shutil.rmtree(self.log_dir)
            params['logging'] = False
        return params

    def get_results_data_params(self):
        params = {}
        reference_filename, target_filename = self.get_results_filenames()
        params['reference'] = self.json_data_dir / reference_filename
        params['target'] = self.json_data_dir / target_filename
        reference_transformed_filename, target_transformed_filename = self.get_transformed_filenames()
        params['reference_relaxed'] = self.json_data_dir / reference_transformed_filename
        params['target_relaxed'] = self.json_data_dir / target_transformed_filename
        return params

    def get_pdb_out_params(self):
        params = {'pdb_out_dir': self.pdb_out_dir}
        for key in self.config['pdb_out_prefix']:
            params.update({key: self.config.get('pdb_out_prefix', key)})
        return params

    def get_results_filenames(self):
        reference_filename = self.config.get('output_files', 'reference_results_file')
        target_filename = self.config.get('output_files', 'target_results_file')
        return reference_filename, target_filename

    def get_transformed_filenames(self):
        reference_filename = self.config.get('output_files', 'transform_reference_results_file')
        target_filename = self.config.get('output_files', 'transform_target_results_file')
        return reference_filename, target_filename

    def get_mutation_params(self):
        params = {}
        try:
            params['repack_radius'] = float(self.config.get('mutation', 'repack_radius'))
        except (configparser.NoSectionError, ValueError):
            params['repack_radius'] = 10.0
        return params

    def get_visualisation_params(self):
        params = {}
        try:
            for key in self.config['visualization']:
                params.update({key: self.config.get('visualization', key)})
            params['out_dir'] = self.out_visualization_dir
        except (configparser.NoSectionError, ValueError):
            params['save'] = False
            shutil.rmtree(self.out_visualization_dir)
        return params

    def clear_tmp_ws(self):
        if os.path.exists(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)
