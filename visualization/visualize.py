import json
import os.path
from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB.Polypeptide import one_to_three
from matplotlib import colors


class VisualizeMng(object):
    def __init__(self,
                 config_params):
        self.reference_data = None
        self.target_data = None
        self.transformed_ref_data = None
        self.transformed_target_data = None
        self.config_params = config_params
        try:
            self.font_size = int(config_params['font_size'])
            self.bar_width = float(config_params['bar_width'])
        except (KeyError, ValueError):
            self.font_size = 10
            self.bar_width = 0.35

        mpl.rcParams.update({'font.size': self.font_size})

    def init_data(self,
                  reference,
                  target,
                  reference_transformed,
                  target_transformed,
                  ):
        self.reference_data = reference
        self.target_data = target
        self.transformed_ref_data = reference_transformed
        self.transformed_target_data = target_transformed

    def load_data(self, data_dir):
        ref_file = os.path.join(data_dir, "reference.json")
        target_file = os.path.join(data_dir, "modified.json")
        transformed_ref_file = os.path.join(data_dir, "reference_modified.json")
        transformed_target_file = os.path.join(data_dir, "reference_relaxed.json")
        self.reference_data = self.load(ref_file)
        self.target_data = self.load(target_file)
        self.transformed_ref_data = self.load(transformed_ref_file)
        self.transformed_target_data = self.load(transformed_target_file)

    @staticmethod
    def load(filepath):
        assert os.path.exists(filepath), f"File not found {filepath}"
        with open(filepath) as fp:
            data = json.load(fp)
        return data

    def plot_residues_energies(self, out_file):
        fig, ax = plt.subplots(nrows=1, ncols=2)
        self.plot_residues_diff(self.reference_data,
                                self.target_data,
                                ax=ax[0],
                                title='reference vs modified',
                                color='royalblue')
        self.plot_residues_diff(self.transformed_ref_data,
                                self.transformed_target_data,
                                ax=ax[1],
                                title='reference (relaxed) vs modified (relaxed)',
                                color='royalblue')
        fig.suptitle('delta G difference')
        plt.savefig(out_file)

    def plot_residues_diff(self,
                           reference_data,
                           target_data,
                           ax,
                           title,
                           color,
                           alpha=None):
        assert self.reference_data is not None
        assert self.target_data is not None
        ref_energies = reference_data['residues_energy']
        ref_x, ref_y = self.__get_total_scores(ref_energies)
        target_energies = target_data['residues_energy']
        target_x, target_y = self.__get_total_scores(target_energies)
        assert len(target_x) == len(ref_x)
        ref_y = np.array(ref_y)
        target_y = np.array(target_y)
        diff = target_y - ref_y
        ref_x = np.array(ref_x, dtype=np.int32)
        ax.set_xlabel('Residues')
        ax.set_ylabel('delta G')
        self.plot_mutation_points_annotation(reference_data, ax, ref_x, diff)
        ax.plot(ref_x, diff, color=color)
        ax.set_title(title)

    @staticmethod
    def plot_mutation_points_annotation(data, ax, ref_x, ref_y):
        mutation_points = data['mutation_points']
        sequence = data['sequence']
        for idx in mutation_points:
            point = int(idx)
            annot_txt = f'Mutation point: {point} {one_to_three(sequence[point - 1])}->{mutation_points[idx]}'
            ax.annotate(annot_txt, xy=(ref_x[point], ref_y[point]),
                        arrowprops=dict(facecolor='r', linewidth=1.0, shrink=0.05))

    def plot_scores(self, out_file):
        sf_fa = []
        sf_ws = []
        sf_weighted = []
        rg = []
        sf_fa_transformed = []
        sf_ws_transformed = []
        sf_weighted_transformed = []
        rg_transformed = []
        titles = ['reference vs modified', 'reference (relaxed) vs modified (relaxed)']
        for rec in (self.reference_data, self.target_data):
            sf_fa.append(rec['sf_fa'])
            sf_ws.append(rec['sf_ws'])
            sf_weighted.append(rec['sf_weighted'])
            # rg.append(rec['radius_gyration'])

        for rec in (self.transformed_ref_data, self.transformed_target_data):
            sf_fa_transformed.append(rec['sf_fa'])
            sf_ws_transformed.append(rec['sf_ws'])
            sf_weighted_transformed.append(rec['sf_weighted'])
            # rg_transformed.append(rec['radius_gyration'])

        plot_rec = (sf_fa, sf_weighted, sf_ws)  # rg)
        plot_rec_transformed = (sf_fa_transformed, sf_weighted_transformed, sf_ws_transformed)  # , rg_transformed)
        fig, ax = plt.subplots(nrows=1, ncols=2)
        for i, item in enumerate((plot_rec, plot_rec_transformed)):
            self.plot_score_hist(item, ax[i], title=titles[i])
        fig.tight_layout()
        plt.savefig(out_file)

    def plot_score_hist(self, data, ax, title):
        labels = ['full atom score', 'ref2015-docking score', 'weighted score']  # , 'radius of gyration']
        ref_values = []
        target_values = []
        for item in data:
            ref, target = item
            ref_values.append(ref)
            target_values.append(target)

        ref_values = np.round(ref_values, 1)
        target_values = np.round(target_values, 1)
        x = np.arange(len(labels))
        width = self.bar_width
        r1 = ax.bar(x - width / 2, ref_values, width, label='reference')
        r2 = ax.bar(x + width / 2, target_values, width, label='modified')
        ax.set_ylabel('Scores')
        ax.set_title(label=title)
        ax.set_xticks(x, labels, rotation='vertical')
        ax.bar_label(r1, padding=-self.font_size * 2, rotation='vertical')
        ax.bar_label(r2, padding=-self.font_size * 2, rotation='vertical')
        ax.legend()

    @staticmethod
    def get_bound_hmap(h_bonds_record, full_size):
        x = []
        y = []
        z = []
        for res_idx in h_bonds_record:
            bound_rec = h_bonds_record.get(res_idx)
            if not bound_rec:
                continue
            acc_idx = bound_rec['acceptor_idx']
            donor_idx = bound_rec['donor_idx']
            score = bound_rec['bound_score']
            x.append(acc_idx)
            y.append(donor_idx)
            z.append(score)
        h_map = np.zeros((full_size, full_size))
        comb = zip(x, y)
        for i, item in enumerate(comb):
            xi, yi = item
            h_map[xi - 1, yi - 1] = 1.0
        return h_map

    @staticmethod
    def __get_total_scores(energies_rec):
        x = []
        y = []
        for res_idx in energies_rec:
            x.append(res_idx)
            y.append(energies_rec[res_idx]['total_score'])
        x = np.array(x, dtype=int)
        y = np.array(y, dtype=float)
        return x, y

    def plot_contacts(self, out_file):
        cm_colors = [['b', 'm'], ['b', 'm']]
        titles = [['reference', 'modified'], ['reference + relaxed', 'modified + relaxed']]
        fig, ax = plt.subplots(nrows=2, ncols=2)
        full_size = len(self.reference_data['sequence'])
        cm_ref = self.get_contact_map(self.reference_data['contact_map'], full_size)
        cm_target = self.get_contact_map(self.target_data['contact_map'], full_size)
        cm_ref_transformed = self.get_contact_map(self.transformed_ref_data['contact_map'], full_size)
        cm_target_transformed = self.get_contact_map(self.transformed_ref_data['contact_map'], full_size)
        cm_rec = [[cm_ref, cm_target], [cm_ref_transformed, cm_target_transformed]]
        for i, cm_item in enumerate(cm_rec):
            for j, cm in enumerate(cm_item):
                self.__plot_hmap_(cm, ax[i][j],
                                  fg_color=cm_colors[i][j],
                                  title=titles[i][j],
                                  axis_labels=('Residues', 'Residues'))
                self.plot_mutation_points_annotation(self.reference_data, ax[i][j],
                                                     range(len(cm.diagonal())), range(len(cm.diagonal())))

        fig.tight_layout()
        fig.suptitle('Contacts map')
        plt.savefig(out_file)

    def plot_h_bonds(self, out_file):
        cm_colors = [['b', 'm'], ['b', 'm']]
        titles = [['reference', 'modified'], ['reference + relaxed', 'modified + relaxed']]
        fig, ax = plt.subplots(nrows=2, ncols=2)
        full_size = len(self.reference_data['sequence'])
        cm_ref = self.get_bound_hmap(self.reference_data['h_bonds'], full_size)
        cm_target = self.get_bound_hmap(self.target_data['h_bonds'], full_size)
        cm_ref_transformed = self.get_bound_hmap(self.transformed_ref_data['h_bonds'], full_size)
        cm_target_transformed = self.get_bound_hmap(self.transformed_ref_data['h_bonds'], full_size)
        cm_rec = [[cm_ref, cm_target], [cm_ref_transformed, cm_target_transformed]]
        for i, cm_item in enumerate(cm_rec):
            for j, cm in enumerate(cm_item):
                self.__plot_hmap_(cm, ax[i][j],
                                  fg_color=cm_colors[i][j],
                                  title=titles[i][j],
                                  axis_labels=('Residues', 'Residues'))

                self.plot_mutation_points_annotation(self.reference_data, ax[i][j],
                                                     range(len(cm.diagonal())), range(len(cm.diagonal())))

        fig.tight_layout()
        fig.suptitle('Hydrogen bonds map')
        plt.savefig(out_file)

    @staticmethod
    def __plot_hmap_(hm,
                     ax,
                     fg_color,
                     bg_color='white',
                     reverse_axis=True,
                     axis_labels=(None, None),
                     title=None):
        color_map = colors.ListedColormap([bg_color] + [fg_color] * 0xf)
        ax.imshow(hm, color_map, interpolation='bicubic')
        if reverse_axis:
            ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_xlabel(axis_labels[0])
        ax.set_ylabel(axis_labels[1])
        ax.set_title(title)

    @staticmethod
    def get_contact_map(contact_record, full_size):
        contact_map = np.zeros((full_size, full_size))
        for res_idx in contact_record:
            contact_list = contact_record[res_idx]
            for res_to in contact_list:
                contact_map[int(res_idx) - 1, int(res_to) - 1] = 1.0

        return contact_map

    def plot_energies_density(self, out_file):
        titles = [['reference', 'modified'], ['reference + relaxed', 'modified + relaxed']]
        fig, ax = plt.subplots(nrows=2, ncols=2)
        density_ref = self.__get_total_scores(self.reference_data['residues_energy'])
        density_target = self.__get_total_scores(self.target_data['residues_energy'])
        density_ref_transformed = self.__get_total_scores(self.transformed_ref_data['residues_energy'])
        density_target_transformed = self.__get_total_scores(self.transformed_ref_data['residues_energy'])
        density_rec = [[density_ref, density_target], [density_ref_transformed, density_target_transformed]]
        for i, density_item in enumerate(density_rec):
            for j, density in enumerate(density_item):
                self.__plot_density(density,
                                    ax[i][j],
                                    axis_labels=('Residues', 'G'),
                                    title=titles[i][j]
                                    )
                self.plot_mutation_points_annotation(self.reference_data, ax[i][j], density[0], density[1])
        plt.autoscale(True)
        fig.suptitle('G-distribution')
        fig.tight_layout()
        plt.savefig(out_file)

    @staticmethod
    def __plot_density(data, ax, axis_labels, title=None):
        normalize = mpl.colors.Normalize(vmin=-1, vmax=1)
        ax.scatter(data[0],
                   data[1],
                   c=data[1],
                   norm=normalize,
                   cmap='jet',
                   marker='o',
                   alpha=0.4,
                   s=8)
        ax.set_xlabel(axis_labels[0])
        ax.set_ylabel(axis_labels[1])
        ax.axhline(y=0.0, color='black', linestyle='--')
        ax.set_title(title)

    def plot(self):
        def __path_resolver(filename):
            return out_dir / f"{filename}.{self.config_params['format']}"

        out_dir = self.config_params['out_dir']
        assert os.path.exists(out_dir), f"Can't access to visualization output {out_dir}"
        out_dir = Path(out_dir)
        out_file = __path_resolver(self.config_params['dg_dif_file'])
        self.plot_residues_energies(out_file)
        out_file = __path_resolver(self.config_params['h_bound_file'])
        self.plot_h_bonds(out_file)
        out_file = __path_resolver(self.config_params['contact_map_file'])
        self.plot_contacts(out_file)
        out_file = __path_resolver(self.config_params['distribution_file'])
        self.plot_energies_density(out_file)
        out_file = __path_resolver(self.config_params['score_file'])
        self.plot_scores(out_file)
        out_file = __path_resolver(self.config_params['ddg_score_file'])
        self.plot_ddg_scores(out_file)

    def plot_ddg_scores(self, out_file):
        sf_fa = []
        sf_ws = []
        sf_weighted = []
        sf_fa_transformed = []
        sf_ws_transformed = []
        sf_weighted_transformed = []
        bar_colors = ['royalblue', 'darkorange', 'g']
        labels = ['full atom', 'ref_2015 docking', 'weighted']

        # ddG = G(unfolded) - G(folded))wt - G(unfolded) - G(folded))mutated

        for rec in (self.reference_data, self.transformed_ref_data):
            sf_fa.append(rec['sf_fa'])
            sf_ws.append(rec['sf_ws'])
            sf_weighted.append(rec['sf_weighted'])

        wt_rec = (sf_fa, sf_ws, sf_ws)
        dG_wt = np.array([g[0] - g[1] for g in wt_rec])

        for rec in (self.target_data, self.transformed_target_data):
            sf_fa_transformed.append(rec['sf_fa'])
            sf_ws_transformed.append(rec['sf_ws'])
            sf_weighted_transformed.append(rec['sf_weighted'])

        folded_rec = (sf_fa_transformed, sf_ws_transformed, sf_ws_transformed)
        dG_folded = np.array([g[0] - g[1] for g in folded_rec])
        ddg = dG_wt - dG_folded
        x = np.arange(0, len(ddg))
        fig, ax = plt.subplots()
        ax.bar(x, ddg, color=bar_colors, width=self.bar_width, label=labels)
        ax.set_xticks(x, labels,  rotation='vertical')
        ax.legend()
        ax.set_title('ddG Score')
        fig.tight_layout()
        plt.savefig(out_file)


test_dir = "/home/dp/Data/mutation/working/6b8eab07-8e83-4267-9ef7-b4ae0892ab2f/results/json_data"

if __name__ == "__main__":
    vis_out = VisualizeMng({})
    vis_out.load_data(test_dir)
    vis_out.plot_ddg_scores('1.jpg')
