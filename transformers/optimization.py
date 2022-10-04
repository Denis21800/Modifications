from pyrosetta.rosetta.protocols.minimization_packing import MinMover
import pyrosetta


class Minimization(object):
    def __init__(self,
                 score_fn):
        mm = pyrosetta.MoveMap()
        mm.set_bb(True)
        mm.set_chi(True)
        self.min_mover = MinMover()
        self.min_mover.movemap(mm)
        self.min_mover.score_function(score_fn)

    def minimize_std(self, pose):
        target_pose = pyrosetta.Pose()
        target_pose.assign(pose)
        self.min_mover.apply(target_pose)
        return target_pose
