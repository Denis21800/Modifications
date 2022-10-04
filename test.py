import sys

from modifiations_tools.mutaion_mngr import MutationMngr, rosetta_configure

TEST_PDB_FILE_1 = '/home/dp/Data/mutation/ms_set/2qug.pdb'
MUTATE_RES_1 = [(155, 'ALY')]

TEST_PDB_FILE_2 = '/home/dp/Data/PDB/1zed.pdb'
MUTATE_RES_2 = [(92, 'SEP')]


if __name__ == "__main__":
    rosetta_configure()
    mngr = MutationMngr()
    task_id, err = mngr.create_task(
        pdb_path=TEST_PDB_FILE_2,
        model_id=0,
        chain_id='A',
        mutate_residues=MUTATE_RES_2
    )
    if task_id == -1:
        print(f"Task error: {err}")
        sys.exit()
    else:
        print(f"Task: {task_id}")
    mngr.mutate(task_id)
