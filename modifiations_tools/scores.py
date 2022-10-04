import json

from pyrosetta import ScoreFunction, create_score_function, get_fa_scorefxn
from pyrosetta.rosetta import core, protocols
from math import sqrt, pow

DEFAULT_CONTACT_CUTOFF = 5.0


class ScoreMetrics(object):
    def __init__(self,
                 pose,
                 score_fn,
                 mutation_points,
                 ):
        self.pose = pose
        self.score_fn = score_fn
        self.results = {'mutation_points': mutation_points}

    def calculate(self):
        self.results['sequence'] = self.pose.sequence()
        fa_score_fxn = get_fa_scorefxn()  # standard full atom ScoreFunction
        ws_patch_score_fxn = create_score_function('ref2015', 'docking')
        score_fxn = ScoreFunction()  # custom
        score_fxn.set_weight(core.scoring.fa_atr, 0.800)  # full-atom attractive score
        score_fxn.set_weight(core.scoring.fa_rep, 0.440)  # full-atom repulsive score
        score_fxn.set_weight(core.scoring.fa_sol, 0.750)  # full-atom solvation score
        score_fxn.set_weight(core.scoring.fa_intra_rep, 0.004)  # f.a. intraresidue rep. score
        score_fxn.set_weight(core.scoring.fa_elec, 0.700)  # full-atom electronic score
        score_fxn.set_weight(core.scoring.pro_close, 1.000)  # proline closure
        score_fxn.set_weight(core.scoring.hbond_sr_bb, 1.170)  # short-range hbonding
        score_fxn.set_weight(core.scoring.hbond_lr_bb, 1.170)  # long-range hbonding
        score_fxn.set_weight(core.scoring.hbond_bb_sc, 1.170)  # backbone-sidechain hbonding
        score_fxn.set_weight(core.scoring.hbond_sc, 1.100)  # sidechain-sidechain hbonding
        score_fxn.set_weight(core.scoring.dslf_fa13, 1.000)  # disulfide full-atom score
        score_fxn.set_weight(core.scoring.rama, 0.200)  # ramachandran score
        score_fxn.set_weight(core.scoring.omega, 0.500)  # omega torsion score
        score_fxn.set_weight(core.scoring.fa_dun, 0.560)  # full atom Dunbrack rotamer score
        score_fxn.set_weight(core.scoring.p_aa_pp, 0.320)
        score_fxn.set_weight(core.scoring.ref, 1.000)  # reference identity score

        self.results['sf_fa'] = fa_score_fxn(self.pose)  # full atom score function
        self.results['sf_ws'] = ws_patch_score_fxn(self.pose)  # full ref2015 docking
        self.results['sf_weighted'] = score_fxn(self.pose)  # custom weighted
        energies = self.pose.energies()
        residue_energies = [energies.residue_total_energy(i)
                            for i in range(1, self.pose.total_residue() + 1)]

        weights = [core.scoring.ScoreType(s)
                   for s in range(1, int(core.scoring.end_of_score_type_enumeration) + 1)
                   if score_fxn.weights()[core.scoring.ScoreType(s)]]

        residue_weighted_energies_matrix = [
            [energies.residue_total_energies(i)[w] * score_fxn.weights()[w]
             for i in range(1, self.pose.total_residue() + 1)]
            for w in weights]

        pose_h_bonds = core.scoring.hbonds.HBondSet()
        core.scoring.hbonds.fill_hbond_set(self.pose, False, pose_h_bonds)
        h_bond_dictionary = {}
        for residue in range(1, self.pose.total_residue() + 1):
            h_bond_record = {}
            for h_bond in range(1, pose_h_bonds.nhbonds() + 1):
                h_bond = pose_h_bonds.hbond(h_bond)
                acceptor_residue = h_bond.acc_res()
                donor_residue = h_bond.don_res()
                if residue == acceptor_residue or residue == donor_residue:
                    h_bond_record['acceptor_idx'] = acceptor_residue
                    h_bond_record['donor_idx'] = donor_residue
                    h_bond_record['acceptor_atom'] = self.pose.residue(acceptor_residue).atom_name(
                        h_bond.acc_atm()).strip()
                    h_bond_record['donor_atom'] = self.pose.residue(donor_residue).atom_name(
                        h_bond.don_hatm()).strip()
                    h_bond_record['bound_score'] = h_bond.energy()

            h_bond_dictionary[residue] = h_bond_record
        self.results['h_bonds'] = h_bond_dictionary

        RadG = ScoreFunction()
        RadG.set_weight(core.scoring.rg, 1)
        pose_radius_gyration = RadG(self.pose)
        self.results['radius_gyration'] = pose_radius_gyration
        residue_energy = {}
        for i in range(1, self.pose.total_residue() + 1):
            residue_detail = {}
            for w in range(len(weights)):
                type_ = core.scoring.name_from_score_type(weights[w])
                residue_detail.update({type_: residue_weighted_energies_matrix[w][i - 1]})

            residue_rec = {'total_score': residue_energies[i - 1], 'weighted': residue_detail}
            residue_energy.update({i: residue_rec})

        self.results['residues_energy'] = residue_energy
        self.results['contact_map'] = self.get_full_contact_map(self.pose)

    def get_score_metrics(self):
        return self.results

    def dump(self, filename):
        with open(filename, 'w') as fp:
            json.dump(self.results, fp)

    def get_full_contact_map(self, pose, cutoff=DEFAULT_CONTACT_CUTOFF):
        contact_map = {}

        for seq_pos in range(1, pose.size() + 1):
            contacts = []
            center_1 = list(pose.residue(seq_pos).nbr_atom_xyz())
            for residue in pose:
                seq_pos_2 = residue.seqpos()
                if seq_pos != seq_pos_2:
                    center_2 = list(residue.nbr_atom_xyz())

                    if self.__calc_distance(center_1, center_2) < cutoff:
                        contacts.append(seq_pos_2)

            if len(contacts) != 0:
                contact_map[seq_pos] = contacts

        return contact_map

    @staticmethod
    def __calc_distance(vec1, vec2):
        x1 = vec1[0]
        y1 = vec1[1]
        z1 = vec1[2]
        x2 = vec2[0]
        y2 = vec2[1]
        z2 = vec2[2]
        dist = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2))

        return dist
